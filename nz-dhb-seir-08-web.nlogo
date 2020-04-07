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

  log-header-file-name
  log-file-name
  full-file-name
  date-time
  model-name

  alert-level-changes
]

to setup
  clear-all

  ask patches [set pcolor sky + 1]

  set-default-shape locales "circle"
  set-default-shape levels "square 3"


  if use-seed? [random-seed seed]

  set R0-levels read-from-string alert-levels-R0
  set flow-levels read-from-string alert-levels-flow
  set trigger-levels read-from-string alert-level-triggers

  ifelse initialise-from-nz-data? [
    initialise-locales-from-string dhbs
    initialise-connections-from-string connectivity
  ]
  [
    initialise-locales-parametrically
    initialise-connections-parametrically
  ]
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

  set model-name "distributed-seir-08"
  set date-time replace date-and-time "." ":"
  let base-file-name (word model-name "-" date-time)
  set log-file-name (word base-file-name ".csv")
  set full-file-name (word log-folder "/" log-file-name)
  set log-header-file-name (word log-folder "/" base-file-name ".header")
;  if log-all-locales? [
;    file-open log-header-file-name
;    file-print log-file-header
;    file-close
;    file-open full-file-name
;    file-print output-locale-header
;    file-close
;  ]
;
  redraw
  reset-ticks
end

to-report replace [s a b]
  let i position a s
  report replace-item i s b
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

;  if log-all-locales? [
;    file-open full-file-name
;    ask locales [
;      file-print output-locale
;    ]
;    file-close
;  ]

  tick
end


to change-alert-levels
  if alert-policy = "static" [ stop ]

  if alert-policy = "local-random" [
    ask locales [
      let alert-level-change one-of (list -1 0 1)
      let new-level alert-level + alert-level-change
      set-transmission-parameters clamp new-level 0 4
      set alert-level-changes alert-level-changes + 1
    ]
    enact-new-levels
    stop
  ]

  if alert-policy = "local" [
    ask locales [
      let recent-new-tests recent-total new-tests
      if recent-new-tests > 0 [
        let local-rate recent-total new-tests-positive / recent-new-tests
        let a get-new-alert-level local-rate
        if a != alert-level [
          set alert-level-changes alert-level-changes + 1
        ]
        ifelse a < alert-level [
          set alert-level alert-level - 1
        ]
        [
          set alert-level a
        ]
      ]
      set-transmission-parameters alert-level
    ]
    enact-new-levels
    stop
  ]

  if alert-policy = "global" [
    let recent-new-tests recent-total all-new-tests
    if recent-new-tests > 0 [
      let global-rate recent-total all-new-tests-positive / recent-new-tests
      let a get-new-alert-level global-rate
      if a != first [alert-level] of locales [
        set alert-level-changes alert-level-changes + 1
      ]
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
to initialise-locales-parametrically
  let pop-mean population / num-locales
  let pop-var pop-sd-multiplier * pop-sd-multiplier * pop-mean * pop-mean

  let alpha (pop-mean * pop-mean / pop-var)
  let lambda (pop-mean / pop-var)

  create-locales num-locales [
    let xy random-xy-with-buffer 0.05
    setxy item 0 xy item 1 xy
    set pop-0 ceiling (random-gamma alpha lambda)
  ]
  let adjustment population / sum [pop-0] of locales
  ask locales [
    set pop-0 round (pop-0 * adjustment)
  ]
  initialise-locales
end

to-report random-xy-with-buffer [buffer]
  let min-x (buffer / 2) * world-width
  let min-y (buffer / 2) * world-height
  let xrange world-width * (1 - buffer)
  let yrange world-height * (1 - buffer)
  report (list (min-x + random-float xrange) (min-y + random-float yrange))
end


to initialise-locales
  let mean-pop-0 mean [pop-0] of locales
  ask locales [
    set size 5 * sqrt (pop-0 / mean-pop-0)
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

;to initialise-locales-from-file [fname]
;  file-open fname
;  let na file-read-line
;  create-locales 20 []
;  foreach sort locales [ x ->
;    ask x [locale-from-line]
;  ]
;  file-close
;  initialise-locales
;end
;
;to locale-from-line
;  let na file-read
;  set label file-read
;  set label-color black
;  repeat 2 [set na file-read]
;  set pop-0 file-read
;  setxy file-read - 119 file-read - 476
;end


to initialise-locales-from-string [s]
  create-locales 20
  (foreach (but-first split-string s "\n") sort locales [ [line x] ->
    let parameters split-string line " "
    ask x [
      set label item 1 parameters
      set label-color black
      set pop-0 read-from-string item 4 parameters
      setxy (read-from-string item 5 parameters - 119) (read-from-string item 6 parameters - 476)
    ]
  ])
  initialise-locales
end


to setup-levels
  ask locales [
    let x nobody
    hatch 1 [
      set size size * sqrt 2
      set breed levels
      set my-locale myself
      set x self
      set label ""
    ]
    set my-alert-indicator x
  ]
end


;; their connections
to initialise-connections-parametrically
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


;to initialise-connections-from-file [fname]
;  file-open fname
;  let na file-read-line
;  while [ not file-at-end? ] [
;    let v1 locale (file-read - 1)
;    let v2 locale (file-read - 1)
;    let weight file-read
;    ask v1 [
;      create-connection-to v2 [
;        initialise-connection weight
;      ]
;    ]
;    repeat 2 [ set na file-read ]
;  ]
;  file-close
;  reweight-connections
;end


to initialise-connections-from-string [s]
  let edges but-first split-string s "\n"
  foreach edges [ edge ->
    let parameters split-string edge " "
    let v1 locale (read-from-string item 0 parameters - 1)
    let v2 locale (read-from-string item 1 parameters - 1)
    let weight read-from-string item 2 parameters
    ask v1 [
      create-connection-to v2 [
        initialise-connection weight
      ]
    ]
  ]
  reweight-connections
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
  set color orange
  set w 1 / d
  set thickness w
end

to reweight-connections
  ask locales [
    let total-w sum [w] of my-out-connections
    let w-correction 1 / total-w
    let thickness-correction 2 / total-w
    ask my-out-connections [
      set w w * w-correction
      set thickness thickness * thickness-correction
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

to-report string-as-list [str]
  report n-values length str [i -> item i str]
end

to-report split-string [str sep]
  let words []
  let this-word ""
  foreach (string-as-list str) [ c ->
    ifelse c = sep
    [ set words sentence words this-word
      set this-word "" ]
    [ set this-word word this-word c ]
  ]
  ifelse this-word = ""
  [ report words ]
  [ report sentence words this-word ]
end

to-report join-list [lst sep]
  report reduce [ [a b] -> (word a sep b) ] lst
end

;to-report output-locale
;  report join-list (list ticks who pop-0 susceptible exposed presymptomatic infected recovered dead tests tests-positive new-exposed new-presymptomatic new-infected new-recovered new-dead alert-level) ","
;end
;
;to-report output-locale-header
;  report "ticks,who,pop-0,susceptible,exposed,presymptomatic,infected,recovered,dead,tests,tests-positive,new-exposed,new-presymptomatic,new-infected,new-recovered,new-dead,alert-level"
;end
;
;to-report log-file-header
;  let parameters (list "name,value")
;  set parameters lput join-list (list "model.name" model-name) "," parameters
;  set parameters lput join-list (list "log.file.name" log-file-name) "," parameters
;  set parameters lput join-list (list "log.all.locales?" log-all-locales?) "," parameters
;  set parameters lput join-list (list "log.folder" log-folder) "," parameters
;
;  set parameters lput join-list (list "use.seed?" use-seed?) "," parameters
;  set parameters lput join-list (list "seed" seed) "," parameters
;
;  set parameters lput join-list (list "initial.infected" initial-infected) "," parameters
;  set parameters lput join-list (list "uniform.by.pop?" uniform-by-pop?) "," parameters
;  set parameters lput join-list (list "init.alert.level" init-alert-level) "," parameters
;  set parameters lput join-list (list "alert.policy" alert-policy) "," parameters
;  set parameters lput join-list (list "start.lifting.quarantine" start-lifting-quarantine) "," parameters
;  set parameters lput join-list (list "time.horizon" time-horizon) "," parameters
;
;  set parameters lput join-list (list "alert.levels.R0" alert-levels-R0) "," parameters
;  set parameters lput join-list (list "alert.levels.flow" alert-levels-flow) "," parameters
;
;  set parameters lput join-list (list "testing.rate.symptomatic" testing-rate-symptomatic) "," parameters
;  set parameters lput join-list (list "testing.rate.presymptomatic" testing-rate-presymptomatic) "," parameters
;  set parameters lput join-list (list "testing.rate.general" testing-rate-general) "," parameters
;
;  set parameters lput join-list (list "population" population) "," parameters
;  set parameters lput join-list (list "num.locales" num-locales) "," parameters
;  set parameters lput join-list (list "pop.sd.multiplier" pop-sd-multiplier) "," parameters
;  set parameters lput join-list (list "flow.rate" flow-rate) "," parameters
;
;  set parameters lput join-list (list "exposed.to.presymptomatic" exposed-to-presymptomatic) "," parameters
;  set parameters lput join-list (list "presymptomatic.to.infected" presymptomatic-to-infected) "," parameters
;  set parameters lput join-list (list "relative.infectiousness.presymptomatic" relative-infectiousness-presymptomatic) "," parameters
;  set parameters lput join-list (list "infected.to.recovered" infected-to-recovered) "," parameters
;  set parameters lput join-list (list "p.hosp" p-hosp) "," parameters
;  set parameters lput join-list (list "p.icu" p-icu) "," parameters
;  set parameters lput join-list (list "icu.cap" icu-cap) "," parameters
;  set parameters lput join-list (list "cfr.0" cfr-0) "," parameters
;  set parameters lput join-list (list "cfr.1" cfr-1) "," parameters
;
;  set parameters lput join-list (list "min.pxcor" min-pxcor) "," parameters
;  set parameters lput join-list (list "max.pxcor" max-pxcor) "," parameters
;  set parameters lput join-list (list "min.pycor" min-pycor) "," parameters
;  set parameters lput join-list (list "max.pycor" max-pycor) "," parameters
;
;  report join-list parameters "\n"
;end
@#$#@#$#@
GRAPHICS-WINDOW
529
10
941
607
-1
-1
4.0
1
16
1
1
1
0
0
0
1
0
100
0
146
1
1
1
ticks
30.0

BUTTON
532
623
606
657
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
623
679
657
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
624
748
658
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
50.0
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
696
669
869
702
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
1079
20
1423
242
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
540
710
713
743
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
672
681
705
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
1429
118
1502
163
dead
total-dead
0
1
11

SWITCH
720
709
867
742
uniform-by-pop?
uniform-by-pop?
0
1
-1000

PLOT
1078
253
1422
487
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
1077
499
1424
736
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
1429
22
1504
67
all-infected
total-infected + total-presymptomatic + total-exposed
0
1
11

MONITOR
1429
70
1502
115
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
876
610
1047
643
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

INPUTBOX
874
648
1059
708
log-folder
test
1
0
String

SWITCH
835
756
1061
790
initialise-from-nz-data?
initialise-from-nz-data?
1
1
-1000

INPUTBOX
6
756
536
849
dhbs
ID name x y pop x10km y10km\n1 Northland 1700348.426 6061149.38 178500 170.0348426 606.114938\n2 Waitemata 1749476.385 5947203.201 622350 174.9476385 594.7203201\n3 Auckland 1757594.65 5920685.756 538840 175.759465 592.0685756\n4 Counties.Manukau 1766973.195 5885442.521 557790 176.6973195 588.5442521\n5 Waikato 1810656.987 5816377.136 417130 181.0656987 581.6377136\n6 Lakes 1878016.101 5754014.304 110040 187.8016101 575.4014304\n7 Bay.of.Plenty 1896923.147 5815743.332 236820 189.6923147 581.5743332\n8 Tairawhiti 2038766.624 5713797.967 48965 203.8766624 571.3797967\n9 Taranaki 1699580.398 5662304.784 119640 169.9580398 566.2304784\n10 Hawkes.Bay 1934136.144 5612568.62 165360 193.4136144 561.256862\n11 Whanganui 1784458.549 5579703.334 64565 178.4458549 557.9703334\n12 MidCentral 1816860.065 5524911.027 178240 181.6860065 552.4911027\n13 Hutt.Valley 1763986.614 5437753.032 149270 176.3986614 543.7753032\n14 Capital.and.Coast 1753070.355 5437211.947 316620 175.3070355 543.7211947\n15 Wairarapa 1817576.511 5457737.439 44845 181.7576511 545.7737439\n16 Nelson.Marlborough 1636138.569 5423667.669 150280 163.6138569 542.3667669\n17 West.Coast 1458495.717 5311071.999 32445 145.8495717 531.1071999\n18 Canterbury 1566860.177 5181157.507 562870 156.6860177 518.1157507\n19 South.Canterbury 1455463.883 5085777.75 60025 145.5463883 508.577775\n20 Southern 1343425.477 4917902.664 328230 134.3425477 491.7902664
1
1
String

INPUTBOX
542
756
828
848
connectivity
ID1 ID2 cost weight weight2\n9 10 0.227069 4.403940 19.394684\n9 14 0.193776 5.160594 26.631733\n9 11 0.093423 10.704056 114.576815\n9 8 0.345958 2.890528 8.355151\n9 4 0.198042 5.049430 25.496742\n9 1 0.335123 2.983978 8.904126\n9 6 0.167174 5.981790 35.781810\n9 7 0.199933 5.001679 25.016795\n9 2 0.236637 4.225886 17.858116\n9 3 0.217664 4.594228 21.106932\n9 5 0.147367 6.785774 46.046732\n9 12 0.131898 7.581596 57.480603\n9 13 0.196719 5.083381 25.840766\n9 15 0.191300 5.227401 27.325723\n9 20 0.702723 1.423036 2.025031\n9 19 0.546829 1.828725 3.344233\n9 18 0.452491 2.209987 4.884042\n9 16 0.309359 3.232490 10.448993\n9 17 0.465441 2.148499 4.616047\n10 9 0.227069 4.403940 19.394684\n10 14 0.174988 5.714684 32.657614\n10 11 0.143124 6.986926 48.817133\n10 8 0.154221 6.484195 42.044779\n10 4 0.229013 4.366559 19.066835\n10 1 0.366094 2.731537 7.461295\n10 6 0.119525 8.366462 69.997687\n10 7 0.164249 6.088325 37.067702\n10 2 0.267608 3.736811 13.963754\n10 3 0.248636 4.021951 16.176089\n10 5 0.171767 5.821835 33.893763\n10 12 0.103927 9.622145 92.585678\n10 13 0.176820 5.655482 31.984481\n10 15 0.134911 7.412309 54.942330\n10 20 0.683935 1.462128 2.137818\n10 19 0.528041 1.893793 3.586451\n10 18 0.433703 2.305725 5.316368\n10 16 0.290571 3.441503 11.843945\n10 17 0.446653 2.238875 5.012560\n14 9 0.193776 5.160594 26.631733\n14 10 0.174988 5.714684 32.657614\n14 11 0.109831 9.104883 82.898897\n14 8 0.326629 3.061574 9.373238\n14 4 0.353972 2.825081 7.981080\n14 1 0.491053 2.036439 4.147086\n14 6 0.246502 4.056755 16.457265\n14 7 0.291226 3.433755 11.790674\n14 2 0.392567 2.547337 6.488927\n14 3 0.373594 2.676699 7.164719\n14 5 0.297856 3.357331 11.271675\n14 12 0.074540 13.415666 179.980082\n14 13 0.012081 82.771665 6851.148445\n14 15 0.053941 18.538685 343.682860\n14 20 0.509581 1.962397 3.851002\n14 19 0.353687 2.827358 7.993953\n14 18 0.259349 3.855804 14.867228\n14 16 0.116217 8.604601 74.039161\n14 17 0.272299 3.672432 13.486755\n11 9 0.093423 10.704056 114.576815\n11 10 0.143124 6.986926 48.817133\n11 14 0.109831 9.104883 82.898897\n11 8 0.294766 3.392522 11.509203\n11 4 0.257996 3.876026 15.023576\n11 1 0.395077 2.531151 6.406726\n11 6 0.164657 6.073228 36.884097\n11 7 0.209381 4.775982 22.810001\n11 2 0.296591 3.371649 11.368015\n11 3 0.277618 3.602066 12.974877\n11 5 0.207321 4.823434 23.265513\n11 12 0.047953 20.853585 434.871996\n11 13 0.112774 8.867253 78.628178\n11 15 0.107355 9.314917 86.767676\n11 20 0.618778 1.616088 2.611742\n11 19 0.462884 2.160367 4.667187\n11 18 0.368546 2.713362 7.362335\n11 16 0.225414 4.436280 19.680584\n11 17 0.381496 2.621257 6.870990\n8 9 0.345958 2.890528 8.355151\n8 10 0.154221 6.484195 42.044779\n8 14 0.326629 3.061574 9.373238\n8 11 0.294766 3.392522 11.509203\n8 4 0.270829 3.692362 13.633541\n8 1 0.407780 2.452302 6.013784\n8 6 0.182956 5.465785 29.874806\n8 7 0.157261 6.358875 40.435291\n8 2 0.309294 3.233172 10.453399\n8 3 0.290321 3.444457 11.864286\n8 5 0.225684 4.430975 19.633541\n8 12 0.255568 3.912846 15.310363\n8 13 0.328461 3.044500 9.268983\n8 15 0.286552 3.489765 12.178459\n8 20 0.835576 1.196779 1.432280\n8 19 0.679682 1.471275 2.164651\n8 18 0.585345 1.708395 2.918615\n8 16 0.442212 2.261358 5.113739\n8 17 0.598294 1.671418 2.793637\n4 9 0.198042 5.049430 25.496742\n4 10 0.229013 4.366559 19.066835\n4 14 0.353972 2.825081 7.981080\n4 11 0.257996 3.876026 15.023576\n4 8 0.270829 3.692362 13.633541\n4 1 0.146393 6.830927 46.661562\n4 6 0.127213 7.860810 61.792338\n4 7 0.113836 8.784570 77.168665\n4 2 0.047907 20.873920 435.720551\n4 3 0.028934 34.561005 1194.463096\n4 5 0.057950 17.256279 297.779174\n4 12 0.287624 3.476764 12.087889\n4 13 0.356916 2.801783 7.849990\n4 15 0.341280 2.930144 8.585745\n4 20 0.862919 1.158857 1.342950\n4 19 0.707025 1.414377 2.000461\n4 18 0.612687 1.632154 2.663925\n4 16 0.469555 2.129676 4.535518\n4 17 0.625637 1.598370 2.554787\n1 9 0.335123 2.983978 8.904126\n1 10 0.366094 2.731537 7.461295\n1 14 0.491053 2.036439 4.147086\n1 11 0.395077 2.531151 6.406726\n1 8 0.407780 2.452302 6.013784\n1 4 0.146393 6.830927 46.661562\n1 6 0.264294 3.783661 14.316091\n1 7 0.250787 3.987450 15.899760\n1 2 0.099803 10.019776 100.395905\n1 3 0.118023 8.472919 71.790359\n1 5 0.195031 5.127394 26.290168\n1 12 0.424705 2.354577 5.544034\n1 13 0.493996 2.024306 4.097815\n1 15 0.478361 2.090471 4.370069\n1 20 1.000000 1.000000 1.000000\n1 19 0.844106 1.184685 1.403479\n1 18 0.749768 1.333745 1.778876\n1 16 0.606636 1.648435 2.717338\n1 17 0.762718 1.311100 1.718983\n6 9 0.167174 5.981790 35.781810\n6 10 0.119525 8.366462 69.997687\n6 14 0.246502 4.056755 16.457265\n6 11 0.164657 6.073228 36.884097\n6 8 0.182956 5.465785 29.874806\n6 4 0.127213 7.860810 61.792338\n6 1 0.264294 3.783661 14.316091\n6 7 0.045598 21.931000 480.968778\n6 2 0.165808 6.031074 36.373859\n6 3 0.146836 6.810337 46.380692\n6 5 0.069967 14.292411 204.273008\n6 12 0.180154 5.550807 30.811463\n6 13 0.249446 4.008888 16.071183\n6 15 0.233810 4.276971 18.292478\n6 20 0.755449 1.323716 1.752223\n6 19 0.599555 1.667902 2.781898\n6 18 0.505218 1.979345 3.917806\n6 16 0.362085 2.761780 7.627430\n6 17 0.518168 1.929878 3.724428\n7 9 0.199933 5.001679 25.016795\n7 10 0.164249 6.088325 37.067702\n7 14 0.291226 3.433755 11.790674\n7 11 0.209381 4.775982 22.810001\n7 8 0.157261 6.358875 40.435291\n7 4 0.113836 8.784570 77.168665\n7 1 0.250787 3.987450 15.899760\n7 6 0.045598 21.931000 480.968778\n7 2 0.152300 6.565968 43.111932\n7 3 0.133328 7.500292 56.254375\n7 5 0.068691 14.558031 211.936253\n7 12 0.224878 4.446857 19.774539\n7 13 0.294170 3.399399 11.555911\n7 15 0.278534 3.590222 12.889695\n7 20 0.800173 1.249729 1.561824\n7 19 0.644279 1.552122 2.409081\n7 18 0.549942 1.818375 3.306487\n7 16 0.406809 2.458155 6.042524\n7 17 0.562891 1.776541 3.156099\n2 9 0.236637 4.225886 17.858116\n2 10 0.267608 3.736811 13.963754\n2 14 0.392567 2.547337 6.488927\n2 11 0.296591 3.371649 11.368015\n2 8 0.309294 3.233172 10.453399\n2 4 0.047907 20.873920 435.720551\n2 1 0.099803 10.019776 100.395905\n2 6 0.165808 6.031074 36.373859\n2 7 0.152300 6.565968 43.111932\n2 3 0.019537 51.185638 2619.969560\n2 5 0.096545 10.357917 107.286450\n2 12 0.326218 3.065432 9.396870\n2 13 0.395510 2.528380 6.392707\n2 15 0.379875 2.632447 6.929776\n2 20 0.901514 1.109246 1.230426\n2 19 0.745620 1.341166 1.798726\n2 18 0.651282 1.535433 2.357555\n2 16 0.508150 1.967924 3.872725\n2 17 0.664232 1.505498 2.266525\n3 9 0.217664 4.594228 21.106932\n3 10 0.248636 4.021951 16.176089\n3 14 0.373594 2.676699 7.164719\n3 11 0.277618 3.602066 12.974877\n3 8 0.290321 3.444457 11.864286\n3 4 0.028934 34.561005 1194.463096\n3 1 0.118023 8.472919 71.790359\n3 6 0.146836 6.810337 46.380692\n3 7 0.133328 7.500292 56.254375\n3 2 0.019537 51.185638 2619.969560\n3 5 0.077572 12.891220 166.183542\n3 12 0.307246 3.254721 10.593207\n3 13 0.376538 2.655776 7.053146\n3 15 0.360902 2.770832 7.677511\n3 20 0.882541 1.133091 1.283896\n3 19 0.726648 1.376183 1.893880\n3 18 0.632310 1.581503 2.501153\n3 16 0.489177 2.044248 4.178951\n3 17 0.645260 1.549764 2.401768\n5 9 0.147367 6.785774 46.046732\n5 10 0.171767 5.821835 33.893763\n5 14 0.297856 3.357331 11.271675\n5 11 0.207321 4.823434 23.265513\n5 8 0.225684 4.430975 19.633541\n5 4 0.057950 17.256279 297.779174\n5 1 0.195031 5.127394 26.290168\n5 6 0.069967 14.292411 204.273008\n5 7 0.068691 14.558031 211.936253\n5 2 0.096545 10.357917 107.286450\n5 3 0.077572 12.891220 166.183542\n5 12 0.231507 4.319521 18.658259\n5 13 0.300799 3.324480 11.052167\n5 15 0.285164 3.506760 12.297363\n5 20 0.806802 1.239461 1.536263\n5 19 0.650909 1.536314 2.360260\n5 18 0.556571 1.796716 3.228190\n5 16 0.413438 2.418739 5.850301\n5 17 0.569521 1.755862 3.083052\n12 9 0.131898 7.581596 57.480603\n12 10 0.103927 9.622145 92.585678\n12 14 0.074540 13.415666 179.980082\n12 11 0.047953 20.853585 434.871996\n12 8 0.255568 3.912846 15.310363\n12 4 0.287624 3.476764 12.087889\n12 1 0.424705 2.354577 5.544034\n12 6 0.180154 5.550807 30.811463\n12 7 0.224878 4.446857 19.774539\n12 2 0.326218 3.065432 9.396870\n12 3 0.307246 3.254721 10.593207\n12 5 0.231507 4.319521 18.658259\n12 13 0.077483 12.906049 166.566106\n12 15 0.065790 15.199795 231.033760\n12 20 0.583487 1.713835 2.937232\n12 19 0.427593 2.338674 5.469394\n12 18 0.333255 3.000705 9.004233\n12 16 0.190123 5.259764 27.665114\n12 17 0.346205 2.888463 8.343220\n13 9 0.196719 5.083381 25.840766\n13 10 0.176820 5.655482 31.984481\n13 14 0.012081 82.771665 6851.148445\n13 11 0.112774 8.867253 78.628178\n13 8 0.328461 3.044500 9.268983\n13 4 0.356916 2.801783 7.849990\n13 1 0.493996 2.024306 4.097815\n13 6 0.249446 4.008888 16.071183\n13 7 0.294170 3.399399 11.555911\n13 2 0.395510 2.528380 6.392707\n13 3 0.376538 2.655776 7.053146\n13 5 0.300799 3.324480 11.052167\n13 12 0.077483 12.906049 166.566106\n13 15 0.043050 23.228547 539.565384\n13 20 0.514490 1.943674 3.777867\n13 19 0.358596 2.788654 7.776593\n13 18 0.264258 3.784180 14.320016\n13 16 0.121126 8.255886 68.159652\n13 17 0.277208 3.607400 13.013335\n15 9 0.191300 5.227401 27.325723\n15 10 0.134911 7.412309 54.942330\n15 14 0.053941 18.538685 343.682860\n15 11 0.107355 9.314917 86.767676\n15 8 0.286552 3.489765 12.178459\n15 4 0.341280 2.930144 8.585745\n15 1 0.478361 2.090471 4.370069\n15 6 0.233810 4.276971 18.292478\n15 7 0.278534 3.590222 12.889695\n15 2 0.379875 2.632447 6.929776\n15 3 0.360902 2.770832 7.677511\n15 5 0.285164 3.506760 12.297363\n15 12 0.065790 15.199795 231.033760\n15 13 0.043050 23.228547 539.565384\n15 20 0.556350 1.797431 3.230759\n15 19 0.400456 2.497155 6.235783\n15 18 0.306118 3.266715 10.671429\n15 16 0.162986 6.135514 37.644530\n15 17 0.319068 3.134130 9.822773\n20 9 0.702723 1.423036 2.025031\n20 10 0.683935 1.462128 2.137818\n20 14 0.509581 1.962397 3.851002\n20 11 0.618778 1.616088 2.611742\n20 8 0.835576 1.196779 1.432280\n20 4 0.862919 1.158857 1.342950\n20 1 1.000000 1.000000 1.000000\n20 6 0.755449 1.323716 1.752223\n20 7 0.800173 1.249729 1.561824\n20 2 0.901514 1.109246 1.230426\n20 3 0.882541 1.133091 1.283896\n20 5 0.806802 1.239461 1.536263\n20 12 0.583487 1.713835 2.937232\n20 13 0.514490 1.943674 3.777867\n20 15 0.556350 1.797431 3.230759\n20 19 0.160045 6.248259 39.040736\n20 18 0.251240 3.980256 15.842441\n20 16 0.478771 2.088682 4.362590\n20 17 0.353693 2.827311 7.993686\n19 9 0.546829 1.828725 3.344233\n19 10 0.528041 1.893793 3.586451\n19 14 0.353687 2.827358 7.993953\n19 11 0.462884 2.160367 4.667187\n19 8 0.679682 1.471275 2.164651\n19 4 0.707025 1.414377 2.000461\n19 1 0.844106 1.184685 1.403479\n19 6 0.599555 1.667902 2.781898\n19 7 0.644279 1.552122 2.409081\n19 2 0.745620 1.341166 1.798726\n19 3 0.726648 1.376183 1.893880\n19 5 0.650909 1.536314 2.360260\n19 12 0.427593 2.338674 5.469394\n19 13 0.358596 2.788654 7.776593\n19 15 0.400456 2.497155 6.235783\n19 20 0.160045 6.248259 39.040736\n19 18 0.095346 10.488083 109.999876\n19 16 0.322877 3.097153 9.592358\n19 17 0.197799 5.055632 25.559411\n18 9 0.452491 2.209987 4.884042\n18 10 0.433703 2.305725 5.316368\n18 14 0.259349 3.855804 14.867228\n18 11 0.368546 2.713362 7.362335\n18 8 0.585345 1.708395 2.918615\n18 4 0.612687 1.632154 2.663925\n18 1 0.749768 1.333745 1.778876\n18 6 0.505218 1.979345 3.917806\n18 7 0.549942 1.818375 3.306487\n18 2 0.651282 1.535433 2.357555\n18 3 0.632310 1.581503 2.501153\n18 5 0.556571 1.796716 3.228190\n18 12 0.333255 3.000705 9.004233\n18 13 0.264258 3.784180 14.320016\n18 15 0.306118 3.266715 10.671429\n18 20 0.251240 3.980256 15.842441\n18 19 0.095346 10.488083 109.999876\n18 16 0.228539 4.375615 19.146003\n18 17 0.149045 6.709402 45.016071\n16 9 0.309359 3.232490 10.448993\n16 10 0.290571 3.441503 11.843945\n16 14 0.116217 8.604601 74.039161\n16 11 0.225414 4.436280 19.680584\n16 8 0.442212 2.261358 5.113739\n16 4 0.469555 2.129676 4.535518\n16 1 0.606636 1.648435 2.717338\n16 6 0.362085 2.761780 7.627430\n16 7 0.406809 2.458155 6.042524\n16 2 0.508150 1.967924 3.872725\n16 3 0.489177 2.044248 4.178951\n16 5 0.413438 2.418739 5.850301\n16 12 0.190123 5.259764 27.665114\n16 13 0.121126 8.255886 68.159652\n16 15 0.162986 6.135514 37.644530\n16 20 0.478771 2.088682 4.362590\n16 19 0.322877 3.097153 9.592358\n16 18 0.228539 4.375615 19.146003\n16 17 0.224747 4.449450 19.797603\n17 9 0.465441 2.148499 4.616047\n17 10 0.446653 2.238875 5.012560\n17 14 0.272299 3.672432 13.486755\n17 11 0.381496 2.621257 6.870990\n17 8 0.598294 1.671418 2.793637\n17 4 0.625637 1.598370 2.554787\n17 1 0.762718 1.311100 1.718983\n17 6 0.518168 1.929878 3.724428\n17 7 0.562891 1.776541 3.156099\n17 2 0.664232 1.505498 2.266525\n17 3 0.645260 1.549764 2.401768\n17 5 0.569521 1.755862 3.083052\n17 12 0.346205 2.888463 8.343220\n17 13 0.277208 3.607400 13.013335\n17 15 0.319068 3.134130 9.822773\n17 20 0.353693 2.827311 7.993686\n17 19 0.197799 5.055632 25.559411\n17 18 0.149045 6.709402 45.016071\n17 16 0.224747 4.449450 19.797603
1
0
String

MONITOR
94
624
188
670
alert-activity
alert-level-changes / count locales / ticks
4
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
  <experiment name="locale-sizes-base-rates-under-static-lockdown-levels" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="730"/>
    <enumeratedValueSet variable="testing-rate-presymptomatic">
      <value value="0.025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="testing-rate-general">
      <value value="5.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="use-seed?">
      <value value="true"/>
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
      <value value="&quot;static&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="init-alert-level" first="0" step="1" last="4"/>
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
