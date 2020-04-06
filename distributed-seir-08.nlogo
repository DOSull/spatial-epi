;;
;; Model for exploring spatial epidemiology of COVID-19
;; from the perspective of regionalised control by 'alert levels'
;; governing local R0
;;
;; One of a series of models, see https://github.com/DOSull/spatial-epi
;;
;; David O'Sullivan
;; david.osullivan@vuw.ac.nz
;;


;; localities where control is administered
breed [locales locale]
locales-own [
  pop-0              ;; intial population
  susceptible        ;; susceptible
  exposed            ;; exposed
  presymptomatic     ;; infected but not symptomatic
  infected           ;; infected and symptomatic
  recovered          ;; recovered
  dead               ;; dead

  untested
  tests
  tests-positive

  new-exposed
  new-presymptomatic
  new-infected
  new-recovered
  new-dead

  new-tests
  new-tests-positive

  alert-level        ;; local alert level which controls...
  my-alert-indicator ;; local level turtle indicator
  my-R0              ;; local R0 and
  my-trans-coeff     ;; local transmission rate
]

;; visualization aid to show alert levels
breed [levels level]
levels-own [
  alert-level
  my-locale
]

;; connections to other locales
directed-link-breed [connections connection]
connections-own [
  w ;; inverse distance weighted _outward_ from each
  my-flow-rate
]


globals [
  n-icu
  cfr-tot

  R0-levels         ;; R0s associated with the alert levels
  flow-levels
  trigger-levels

  mean-R0           ;; pop weighted mean R0
  mean-trans-coeff  ;; pop weighted mean trans coeff

  total-exposed            ;; exposed
  total-presymptomatic     ;; infected but not symptomatic
  total-infected           ;; infected and symptomatic
  total-recovered          ;; recovered
  total-dead               ;; dead

  all-tests
  all-tests-positive

  total-new-exposed
  total-new-presymptomatic
  total-new-infected
  total-new-recovered
  total-new-dead

  all-new-tests
  all-new-tests-positive

  log-folder
  log-header-file-name
  log-file-name
  date-time
  model-name
]

to setup
  clear-all

  ask patches [set pcolor blue - 3]

  if use-seed? [random-seed seed]

  set R0-levels read-from-string alert-levels-R0
  set flow-levels read-from-string alert-levels-flow
  set trigger-levels read-from-string alert-level-triggers

  setup-locales
  connect-network
  setup-levels

  ;; initial exposures
  ifelse uniform-by-pop? [
    uniform-expose-by-population initial-infected
  ]
  [
    repeat initial-infected [
      ask one-of locales [
        initial-infect
      ]
    ]
  ]
  set all-new-tests []
  set all-new-tests-positive []
  update-global-parameters

  set model-name "distributed-seir-08-"
  set date-time date-and-time
  set log-folder "distributed-seir-results"
  set log-file-name (word log-folder "/" model-name date-time ".csv")
  set log-header-file-name (word log-folder "/" model-name date-time ".header")
  if log-all-locales? [
    file-open log-header-file-name
    file-print log-file-header
    file-close
    file-open log-file-name
    file-print output-locale-header
    file-close
  ]

  redraw
  reset-ticks
end

;; susceptible weighted infection so
;; that higher population locales are more
;; likely to receive exposures
to uniform-expose-by-population [n]
  repeat n [
    ;; build cumulative total of susceptibles by locale
    ;; ordered from highest to lowest
    let susc reverse sort [susceptible] of locales
    let cum-susc cumulative-sum susc
    let total-susc last cum-susc

    ;; order locales similarly
    let ordered-locales reverse sort-on [susceptible] locales

    ;; pick a random number and use it to pick the corresponding locale
    let i random total-susc
    let idx length filter [ x -> x < i ] cum-susc
    ask item idx ordered-locales [
      initial-infect
    ]
  ]
end


to initial-infect
  set susceptible susceptible - 1

  let choose one-of (list 1 2 3)
  if choose = 1 [
    set exposed exposed + 1
    stop
  ]
  if choose = 2 [
    set presymptomatic presymptomatic + 1
    stop
  ]
  set infected infected + 1
end

to-report cumulative-sum [lst]
  let starter sublist lst 0 1
  report reduce [ [a b] -> lput (last a + b) a] (fput starter but-first lst)
end

to update-global-parameters
  set mean-R0 sum [my-R0 * susceptible] of locales / sum [susceptible] of locales
  set mean-trans-coeff sum [my-trans-coeff * susceptible] of locales / sum [susceptible] of locales

  set total-exposed sum [exposed] of locales
  set total-presymptomatic sum [presymptomatic] of locales
  set total-infected sum [infected] of locales
  set total-recovered sum [recovered] of locales
  set total-dead sum [dead] of locales

  set total-new-exposed sum [new-exposed] of locales
  set total-new-presymptomatic sum [new-presymptomatic] of locales
  set total-new-infected sum [new-infected] of locales
  set total-new-recovered sum [new-recovered] of locales
  set total-new-dead sum [new-dead] of locales

  set all-tests sum [tests] of locales
  set all-tests-positive sum [tests-positive] of locales
  set all-new-tests fput (sum [first new-tests] of locales) all-new-tests
  set all-new-tests-positive fput (sum [first new-tests-positive] of locales) all-new-tests-positive

  set n-icu total-infected * p-icu
  set cfr-tot cfr-0
  if n-icu > icu-cap [
    set cfr-tot (cfr-0 * icu-cap + (n-icu - icu-cap) * cfr-1) / n-icu
  ]
end


;; -------------------------
;; Main
;; -------------------------
to go
  if ((total-infected + total-presymptomatic + total-exposed) = 0 and ticks > 0) [
    update-global-parameters
    update-plots
    tick
    stop
  ]
  update-testing-results
  if ticks >= start-lifting-quarantine and (ticks - start-lifting-quarantine) mod time-horizon = 0 [
    change-alert-levels
  ]
  update-global-parameters

  ask locales [
    spread
  ]
  redraw

  if log-all-locales? [
    file-open log-file-name
    ask locales [
      file-print output-locale
    ]
    file-close
  ]

  tick
end


to change-alert-levels
  if alert-policy = "static" [ stop ]

  if alert-policy = "local-random" [
    ask locales [
      let alert-level-change one-of (list -1 0 1)
      let new-level alert-level + alert-level-change
      set-transmission-parameters clamp new-level 0 4
    ]
    enact-new-levels
    stop
  ]

  if alert-policy = "local" [
    ask locales [
      let recent-new-tests recent-total new-tests
      ifelse recent-new-tests > 0 [
        let local-rate recent-total new-tests-positive / recent-new-tests
        let a get-new-alert-level local-rate
        ifelse a < alert-level [
          set alert-level alert-level - 1
        ]
        [
          set alert-level a
        ]
      ]
      [
        set alert-level 0
      ]
      set-transmission-parameters alert-level
    ]
    enact-new-levels
    stop
  ]

  if alert-policy = "global" [
    let recent-new-tests recent-total all-new-tests
    ifelse recent-new-tests > 0 [
      let global-rate recent-total all-new-tests-positive / recent-new-tests
      let a get-new-alert-level global-rate
      ask locales [
        ifelse a < alert-level [
          set alert-level alert-level - 1
        ]
        [
          set alert-level a
        ]
        set-transmission-parameters alert-level
      ]
    ]
    [
      ask locales [
        set alert-level 0
        set-transmission-parameters alert-level
      ]
    ]
    enact-new-levels
    stop
  ]
end

to-report recent-total [lst]
  report sum sublist lst 0 time-horizon
end

to-report get-new-alert-level [r]
  report length filter [x -> r > x] trigger-levels
end

to-report clamp [x mn mx]
  report max (list mn min (list mx x))
end

to enact-new-levels
  ask levels [
    set alert-level [alert-level] of my-locale
    draw-level
  ]
  ask connections [
    set my-flow-rate min [item alert-level flow-levels] of both-ends
  ]
end


to set-transmission-parameters [a]
  set alert-level a
  set my-R0 item alert-level R0-levels
  set my-trans-coeff get-transmission-coeff my-R0
end

to-report get-transmission-coeff [R]
  report R / (relative-infectiousness-presymptomatic / presymptomatic-to-infected + 1 / infected-to-recovered)
end



to spread
  ;; calculate all the flows
  set new-exposed random-binomial susceptible (my-trans-coeff * (relative-infectiousness-presymptomatic * get-effective-presymptomatic + get-effective-infected) / (pop-0 - dead))
  set new-presymptomatic random-binomial exposed exposed-to-presymptomatic
  set new-infected random-binomial presymptomatic presymptomatic-to-infected

  let no-longer-infected random-binomial infected infected-to-recovered
  set new-recovered random-binomial no-longer-infected (1 - cfr-tot)
  set new-dead no-longer-infected - new-recovered

  ;; update all the stocks
  expose new-exposed
  presym new-presymptomatic
  infect new-infected
  recover new-recovered
  kill new-dead
end

to update-testing-results
  ask locales [
    set new-tests-positive fput (random-binomial infected testing-rate-symptomatic) new-tests-positive

    let new-tests-pre random-binomial (presymptomatic + exposed) testing-rate-presymptomatic
    let new-tests-negative random-binomial susceptible testing-rate-general
    set new-tests fput (first new-tests-positive + new-tests-pre + new-tests-negative) new-tests

    set tests tests + first new-tests
    set tests-positive tests-positive + first new-tests-positive
  ]
end

to expose [n]
  ;;if n != 0 [show word "expose " n]
  set susceptible susceptible - n
  set exposed exposed + n
end

to presym [n]
  ;;if n != 0 [show word "presym " n]
  set exposed exposed - n
  set presymptomatic presymptomatic + n
end

to infect [n]
  ;;if n != 0 [show word "infect " n]
  set presymptomatic presymptomatic - n
  set infected infected + n
end

to recover [n]
  ;;if n != 0 [show word "recover " n]
  set infected infected - n
  set recovered recovered + n
end

to kill [n]
  ;;if n != 0 [show word "kill " n]
  set infected infected - n
  set pop-0 pop-0 - n
  set dead dead + n
end

to-report get-effective-infected
  report infected + flow-rate * sum [my-flow-rate * w * [infected] of other-end] of my-in-connections
end

to-report get-effective-presymptomatic
  report presymptomatic + flow-rate * sum [my-flow-rate * w * [presymptomatic] of other-end] of my-in-connections
end

;; --------------------------------------------------------------
;; NOTE
;; using random-poisson approximation for efficiency when n large
;; --------------------------------------------------------------
to-report random-binomial [n p]
  if n > 100 [report random-poisson (n * p)]
  report length filter [x -> x < p] (n-values n [x -> random-float 1])
end


;; --------------------------------
;; initialisation stuff
;; --------------------------------
;; locales
to setup-locales
  set-default-shape locales "circle"

  let pop-mean population / num-locales
  let pop-var pop-sd-multiplier * pop-sd-multiplier * pop-mean * pop-mean

  let alpha (pop-mean * pop-mean / pop-var)
  let lambda (pop-mean / pop-var)

  create-locales num-locales [
    setxy (random-xcor * 0.95) (random-ycor * 0.95)
    set pop-0 ceiling (random-gamma alpha lambda)
  ]
  let adjustment population / sum [pop-0] of locales
  ask locales [
    set pop-0 round (pop-0 * adjustment)
    set size sqrt (pop-0 / pop-mean)
    set susceptible pop-0
    set exposed 0
    set presymptomatic 0
    set infected 0
    set recovered 0
    set dead 0

    set untested 0
    set tests 0
    set tests-positive 0

    set new-exposed 0
    set new-presymptomatic 0
    set new-infected 0
    set new-recovered 0
    set new-dead 0

    set new-tests (list 0)
    set new-tests-positive (list 0)

    set-transmission-parameters init-alert-level
  ]
end

to setup-levels
  set-default-shape levels "square 3"

  ask locales [
    let x nobody
    hatch 1 [
      set size size * sqrt 2
      set breed levels
      set my-locale myself
      set x self
    ]
    set my-alert-indicator x
  ]
end


;; their connections
to connect-network
  ask locales [
    create-pair-of-connections self nearest-non-neighbour
  ]
  let num-links (6 * count locales)
  while [count connections < num-links] [
    ask one-of locales [
      create-pair-of-connections self nearest-non-neighbour
    ]
  ]
  reweight-connections
end

to-report nearest-non-neighbour
  report first sort-on [distance myself] (other locales with [not connection-neighbor? myself])
end

;; creates directed links between a and b in both directions
to create-pair-of-connections [a b]
  ask a [
    let d distance b
    create-connection-to b [initialise-connection d]
  ]
  ask b [
    let d distance a
    create-connection-to a [initialise-connection d]
  ]
end

to initialise-connection [d]
  set color grey
  set w 1 / d
end

to reweight-connections
  ask locales [
    let total-w sum [w] of my-out-connections
    ask my-out-connections [
      set w w / total-w
      set thickness w / 2
    ]
  ]
end

;; ----------------------------
;; Drawing stuff
;; ----------------------------
to redraw
  ask locales [
    draw-locale
  ]
  ask levels [
    draw-level
  ]
end

to draw-locale
  set color scale-color red dead (cfr-0 * pop-0) 0
end

to draw-level
  set alert-level [alert-level] of my-locale
  set color item alert-level [grey lime yellow orange red]
end


to-report join-list [lst sep]
  report reduce [ [a b] -> (word a sep b) ] lst
end

to-report output-locale
  report join-list (list ticks who pop-0 susceptible exposed presymptomatic infected recovered dead tests tests-positive new-exposed new-presymptomatic new-infected new-recovered new-dead alert-level) ","
end

to-report output-locale-header
  report "ticks,who,pop-0,susceptible,exposed,presymptomatic,infected,recovered,dead,tests,tests-positive,new-exposed,new-presymptomatic,new-infected,new-recovered,new-dead,alert-level"
end

to-report log-file-header
  let parameters (list model-name date-time)
  set parameters lput join-list (list "min-pxcor" min-pxcor) "," parameters
  set parameters lput join-list (list "max-pxcor" max-pxcor) "," parameters
  set parameters lput join-list (list "min-pycor" min-pycor) "," parameters
  set parameters lput join-list (list "max-pycor" max-pycor) "," parameters
  set parameters lput join-list (list "use-seed?" use-seed?) "," parameters
  set parameters lput join-list (list "seed" seed) "," parameters

  set parameters lput join-list (list "initial-infected" initial-infected) "," parameters
  set parameters lput join-list (list "uniform-by-pop?" uniform-by-pop?) "," parameters
  set parameters lput join-list (list "log-all-locales?" log-all-locales?) "," parameters

  set parameters lput join-list (list "init-alert-level" init-alert-level) "," parameters
  set parameters lput join-list (list "alert-policy" alert-policy) "," parameters
  set parameters lput join-list (list "start-lifting-quarantine" start-lifting-quarantine) "," parameters
  set parameters lput join-list (list "time-horizon" time-horizon) "," parameters

  set parameters lput join-list (list "alert-levels-R0" alert-levels-R0) "," parameters
  set parameters lput join-list (list "alert-levels-flow" alert-levels-flow) "," parameters

  set parameters lput join-list (list "testing-rate-symptomatic" testing-rate-symptomatic) "," parameters
  set parameters lput join-list (list "testing-rate-presymptomatic" testing-rate-presymptomatic) "," parameters
  set parameters lput join-list (list "testing-rate-general" testing-rate-general) "," parameters

  set parameters lput join-list (list "population" population) "," parameters
  set parameters lput join-list (list "num-locales" num-locales) "," parameters
  set parameters lput join-list (list "pop-sd-multiplier" pop-sd-multiplier) "," parameters
  set parameters lput join-list (list "flow-rate" flow-rate) "," parameters

  set parameters lput join-list (list "exposed-to-presymptomatic" exposed-to-presymptomatic) "," parameters
  set parameters lput join-list (list "presymptomatic-to-infected" presymptomatic-to-infected) "," parameters
  set parameters lput join-list (list "relative-infectiousness-presymptomatic" relative-infectiousness-presymptomatic) "," parameters
  set parameters lput join-list (list "infected-to-recovered" infected-to-recovered) "," parameters
  set parameters lput join-list (list "p-hosp" p-hosp) "," parameters
  set parameters lput join-list (list "p-icu" p-icu) "," parameters
  set parameters lput join-list (list "icu-cap" icu-cap) "," parameters
  set parameters lput join-list (list "cfr-0" cfr-0) "," parameters
  set parameters lput join-list (list "cfr-1" cfr-1) "," parameters

  report join-list parameters "\n"
end
@#$#@#$#@
GRAPHICS-WINDOW
529
10
1095
577
-1
-1
18.0
1
24
1
1
1
0
0
0
1
-15
15
-15
15
1
1
1
ticks
30.0

BUTTON
532
587
606
621
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
614
587
679
621
step
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
684
588
748
622
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
348
63
522
96
num-locales
num-locales
20
200
100.0
10
1
NIL
HORIZONTAL

SLIDER
12
155
195
188
exposed-to-presymptomatic
exposed-to-presymptomatic
0
1
0.25
0.01
1
NIL
HORIZONTAL

SLIDER
13
196
194
229
presymptomatic-to-infected
presymptomatic-to-infected
0
1
1.0
0.1
1
NIL
HORIZONTAL

SLIDER
12
236
195
269
infected-to-recovered
infected-to-recovered
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
12
277
317
310
relative-infectiousness-presymptomatic
relative-infectiousness-presymptomatic
0
1
0.15
0.01
1
NIL
HORIZONTAL

MONITOR
86
26
213
71
mean-trans-coeff
mean-trans-coeff
5
1
11

SLIDER
14
483
239
516
testing-rate-symptomatic
testing-rate-symptomatic
0
1
0.25
0.01
1
NIL
HORIZONTAL

SLIDER
13
380
153
413
cfr-0
cfr-0
0
0.1
0.01
0.001
1
NIL
HORIZONTAL

SLIDER
159
380
306
413
cfr-1
cfr-1
0
0.2
0.02
0.001
1
NIL
HORIZONTAL

SLIDER
13
338
150
371
p-hosp
p-hosp
0
0.5
0.05
0.01
1
NIL
HORIZONTAL

SLIDER
11
419
153
452
p-icu
p-icu
0
p-hosp
0.0125
0.0001
1
NIL
HORIZONTAL

SLIDER
158
419
308
452
icu-cap
icu-cap
100
600
500.0
10
1
beds
HORIZONTAL

SLIDER
695
633
868
666
initial-infected
initial-infected
0
5000
2500.0
10
1
NIL
HORIZONTAL

PLOT
1107
15
1451
237
totals
days
log (people + 1)
0.0
6.0
0.0
6.0
true
true
"" ""
PENS
"exposed" 1.0 0 -8431303 true "" "plot log (total-exposed + 1) 10"
"presymp" 1.0 0 -955883 true "" "plot log (total-presymptomatic + 1) 10"
"infected" 1.0 0 -2674135 true "" "plot log (total-infected + 1) 10"
"recovered" 1.0 0 -13840069 true "" "plot log (total-recovered + 1) 10"
"dead" 1.0 0 -16777216 true "" "plot log (total-dead + 1) 10"

SLIDER
539
674
712
707
seed
seed
0
100
30.0
1
1
NIL
HORIZONTAL

SLIDER
344
217
517
250
flow-rate
flow-rate
0
1
1.0
0.1
1
NIL
HORIZONTAL

SWITCH
542
636
681
669
use-seed?
use-seed?
1
1
-1000

SLIDER
348
26
521
59
population
population
100000
10000000
5000000.0
100000
1
NIL
HORIZONTAL

TEXTBOX
353
7
437
25
Population\n
12
0.0
1

TEXTBOX
7
6
91
24
Pandemic
12
0.0
1

TEXTBOX
10
297
94
315
Mortality
12
0.0
1

TEXTBOX
345
197
429
215
Connectivity
12
0.0
1

TEXTBOX
12
461
143
491
Control and testing
12
0.0
1

SLIDER
348
101
520
134
pop-sd-multiplier
pop-sd-multiplier
0.01
1.2
0.45
0.01
1
NIL
HORIZONTAL

MONITOR
263
24
343
69
total-pop
sum [pop-0] of locales
0
1
11

MONITOR
437
145
517
190
max-pop
max [pop-0] of locales
0
1
11

MONITOR
350
146
431
191
min-pop
min [pop-0] of locales
0
1
11

INPUTBOX
132
86
321
146
alert-levels-R0
[2.5 2.1 1.6 1.1 0.6]
1
0
String

MONITOR
7
26
79
71
mean-R0
mean-R0
3
1
11

SLIDER
348
287
521
320
init-alert-level
init-alert-level
0
4
4.0
1
1
NIL
HORIZONTAL

MONITOR
1457
113
1530
158
dead
total-dead
0
1
11

SWITCH
720
673
867
706
uniform-by-pop?
uniform-by-pop?
0
1
-1000

PLOT
1106
248
1450
482
new-cases
days
number of people
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"exposed" 1.0 0 -6459832 true "" "plot total-new-exposed"
"presymp" 1.0 0 -955883 true "" "plot total-new-presymptomatic"
"infected" 1.0 0 -2674135 true "" "plot total-new-infected"
"recovered" 1.0 0 -13840069 true "" "plot total-new-recovered"
"dead" 1.0 0 -16777216 true "" "plot total-new-dead"

SLIDER
13
520
237
553
testing-rate-presymptomatic
testing-rate-presymptomatic
0
0.1
0.025
0.001
1
NIL
HORIZONTAL

SLIDER
12
559
236
592
testing-rate-general
testing-rate-general
0
0.002
5.0E-4
1
1
NIL
HORIZONTAL

PLOT
1105
494
1452
731
testing
days
tests
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"tests-pos" 1.0 0 -7500403 true "" "plot all-tests-positive"
"new-tests" 1.0 0 -2674135 true "" "plot first all-new-tests"
"new-tests-pos" 1.0 0 -955883 true "" "plot first all-new-tests-positive"

INPUTBOX
367
519
523
579
alert-levels-flow
[1.0 0.5 0.25 0.1 0.05]
1
0
String

CHOOSER
367
327
520
372
alert-policy
alert-policy
"global" "local" "local-random" "static"
1

MONITOR
1457
17
1532
62
all-infected
total-infected + total-presymptomatic + total-exposed
0
1
11

MONITOR
1457
65
1530
110
recovered
total-recovered
0
1
11

INPUTBOX
316
377
523
437
alert-level-triggers
[0.0005 0.001 0.0025 0.005 1]
1
0
String

SLIDER
349
479
522
512
time-horizon
time-horizon
1
28
7.0
1
1
days
HORIZONTAL

SWITCH
540
714
711
747
log-all-locales?
log-all-locales?
1
1
-1000

SLIDER
322
440
524
473
start-lifting-quarantine
start-lifting-quarantine
0
56
28.0
7
1
days
HORIZONTAL

TEXTBOX
350
269
429
287
Alert levels
12
0.0
1

MONITOR
199
623
257
668
pop-lev-0
sum [pop-0] of locales with [alert-level = 0]
0
1
11

MONITOR
263
623
321
668
pop-lev-1
sum [pop-0] of locales with [alert-level = 1]
0
1
11

MONITOR
325
623
383
668
pop-lev-2
sum [pop-0] of locales with [alert-level = 2]
0
1
11

MONITOR
388
623
447
668
pop-lev-3
sum [pop-0] of locales with [alert-level = 3]
0
1
11

MONITOR
453
623
512
668
pop-lev-4
sum [pop-0] of locales with [alert-level = 4]
0
1
11

MONITOR
199
674
257
719
n-lev-0
count locales with [alert-level = 0]
0
1
11

MONITOR
263
674
321
719
n-lev-1
count locales with [alert-level = 1]
0
1
11

MONITOR
325
674
383
719
n-lev-2
count locales with [alert-level = 2]
0
1
11

MONITOR
388
674
447
719
n-lev-3
count locales with [alert-level = 3]
0
1
11

MONITOR
453
674
512
719
n-lev-4
count locales with [alert-level = 4]
0
1
11

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

square 3
false
0
Rectangle -7500403 false true 15 15 285 285
Rectangle -7500403 false true 0 0 300 300
Rectangle -7500403 false true 7 9 292 293

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="compare-locale-sizes-static" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="730"/>
    <metric>total-infected</metric>
    <metric>total-recovered</metric>
    <metric>total-dead</metric>
    <enumeratedValueSet variable="alert-policy">
      <value value="&quot;static&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="testing-rate-presymptomatic">
      <value value="0.025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="testing-rate-general">
      <value value="5.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-alert-level">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="seed" first="1" step="1" last="30"/>
    <enumeratedValueSet variable="relative-infectiousness-presymptomatic">
      <value value="0.15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-icu">
      <value value="0.0125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pop-sd-multiplier">
      <value value="0.45"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="population">
      <value value="5000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-hosp">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="uniform-by-pop?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-exposed">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infected-to-recovered">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cfr-0">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exposed-to-presymptomatic">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-horizon">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alert-levels-R0">
      <value value="&quot;[2.5 2.1 1.6 1.1 0.6]&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-locales">
      <value value="20"/>
      <value value="50"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cfr-1">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alert-levels-flow">
      <value value="&quot;[1.0 0.5 0.25 0.1 0.05]&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="presymptomatic-to-infected">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alert-level-triggers">
      <value value="&quot;[0.0005 0.001 0.0025 0.005 1]&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="icu-cap">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="testing-rate-symptomatic">
      <value value="0.25"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="compare-locale-sizes-local" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="730"/>
    <metric>total-infected</metric>
    <metric>total-recovered</metric>
    <metric>total-dead</metric>
    <enumeratedValueSet variable="alert-policy">
      <value value="&quot;local&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="testing-rate-presymptomatic">
      <value value="0.025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="testing-rate-general">
      <value value="5.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-alert-level">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="seed" first="1" step="1" last="30"/>
    <enumeratedValueSet variable="relative-infectiousness-presymptomatic">
      <value value="0.15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-icu">
      <value value="0.0125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pop-sd-multiplier">
      <value value="0.45"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="population">
      <value value="5000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-hosp">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="uniform-by-pop?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-exposed">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infected-to-recovered">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cfr-0">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exposed-to-presymptomatic">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-horizon">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alert-levels-R0">
      <value value="&quot;[2.5 2.1 1.6 1.1 0.6]&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-locales">
      <value value="20"/>
      <value value="50"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cfr-1">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alert-levels-flow">
      <value value="&quot;[1.0 0.5 0.25 0.1 0.05]&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="presymptomatic-to-infected">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alert-level-triggers">
      <value value="&quot;[0.0005 0.001 0.0025 0.005 1]&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="icu-cap">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="testing-rate-symptomatic">
      <value value="0.25"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="locale-sizes-vs-lockdown-counts" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="730"/>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="testing-rate-presymptomatic">
      <value value="0.025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="testing-rate-general">
      <value value="5.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use-seed?">
      <value value="false"/>
    </enumeratedValueSet>
    <steppedValueSet variable="seed" first="1" step="1" last="30"/>
    <enumeratedValueSet variable="population">
      <value value="5000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="uniform-by-pop?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infected-to-recovered">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="flow-rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cfr-0">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-horizon">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="exposed-to-presymptomatic">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cfr-1">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alert-levels-flow">
      <value value="&quot;[1.0 0.5 0.25 0.1 0.05]&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="icu-cap">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-lifting-quarantine">
      <value value="28"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="log-all-locales?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-infected">
      <value value="2500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alert-policy">
      <value value="&quot;local&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="init-alert-level">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="relative-infectiousness-presymptomatic">
      <value value="0.15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-icu">
      <value value="0.0125"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pop-sd-multiplier">
      <value value="0.45"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-hosp">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alert-levels-R0">
      <value value="&quot;[2.5 2.1 1.6 1.1 0.6]&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="num-locales">
      <value value="20"/>
      <value value="50"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="presymptomatic-to-infected">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="testing-rate-symptomatic">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alert-level-triggers">
      <value value="&quot;[0.0005 0.001 0.0025 0.005 1]&quot;"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
