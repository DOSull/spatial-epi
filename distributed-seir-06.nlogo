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

  my-flow-rate       ;; not yet implemented: control over inter-regional movement
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
]


globals [
  n-icu
  cfr-tot

  R0-levels         ;; R0s associated with the alert levels
  mean-R0           ;; pop weighted mean R0
  mean-trans-coeff  ;; pop weighted mean trans coeff

  total-exposed            ;; exposed
  total-presymptomatic     ;; infected but not symptomatic
  total-infected           ;; infected and symptomatic
  total-recovered          ;; recovered
  total-dead               ;; dead

  total-tests
  total-tests-positive

  total-new-exposed
  total-new-presymptomatic
  total-new-infected
  total-new-recovered
  total-new-dead

  total-new-tests
  total-new-tests-positive
]

to setup
  clear-all

  ask patches [set pcolor blue - 3]

  if use-seed? [random-seed seed]

  set R0-levels read-from-string alert-levels-R0

  setup-locales
  connect-network
  setup-levels

  ;; initial exposures
  ifelse uniform-by-pop? [
    uniform-expose-by-population initial-exposed
  ]
  [
    repeat initial-exposed [
      ask one-of locales [
        expose 1
      ]
    ]
  ]
  update-global-parameters

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
      expose 1
    ]
  ]
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

  set total-tests sum [tests] of locales
  set total-new-tests sum [new-tests] of locales
  set total-tests-positive sum [tests-positive] of locales
  set total-new-tests-positive sum [new-tests-positive] of locales

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
  if (sum [infected + presymptomatic] of locales = 0 and ticks > 0) [
    stop
  ]
  ;; change-alert-levels
  update-global-parameters

  ask locales [
    spread
  ]
  redraw

  tick
end


to change-alert-levels
  ask locales [
    let alert-level-change one-of (list -1 0 1)
    let new-level alert-level + alert-level-change
    set alert-level clamp new-level 0 4
    ask my-alert-indicator [
      set alert-level [alert-level] of myself
      draw-level
    ]
    set-transmission-parameters alert-level
  ]
end

to-report clamp [x mn mx]
  report max (list mn min (list mx x))
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

  update-testing-results

  ;; update all the stocks
  expose new-exposed
  presym new-presymptomatic
  infect new-infected
  recover new-recovered
  kill new-dead
end

to update-testing-results
  set new-tests-positive random-binomial infected testing-rate-symptomatic
  let new-tests-pre random-binomial (presymptomatic + exposed) testing-rate-presymptomatic
  let new-tests-not random-binomial susceptible testing-rate-general
  set new-tests new-tests-positive + new-tests-pre + new-tests-not

  set tests tests + new-tests
  set tests-positive tests-positive + new-tests-positive
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
  report infected + flow-rate * sum [w * [infected] of other-end] of my-in-connections
end

to-report get-effective-presymptomatic
  report presymptomatic + flow-rate * sum [w * [presymptomatic] of other-end] of my-in-connections
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
    set size sqrt (pop-0 / pop-mean)
    set susceptible pop-0
    set exposed 0
    set presymptomatic 0
    set infected 0
    set recovered 0
    set dead 0

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
  ; make the network look a little prettier
  repeat 5 [
    layout-spring locales connections 0.3 (world-width / (sqrt count locales)) 1
  ]
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
  set color scale-color red dead (cfr-1 * pop-0) 0
end

to draw-level
  set alert-level [alert-level] of my-locale
  set color item alert-level [grey lime yellow orange red]
end

@#$#@#$#@
GRAPHICS-WINDOW
529
10
1096
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
334
66
508
99
num-locales
num-locales
10
200
100.0
10
1
NIL
HORIZONTAL

SLIDER
18
82
190
115
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
19
123
188
156
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
18
163
188
196
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
18
204
323
237
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
98
26
225
72
mean-trans-coeff
mean-trans-coeff
5
1
11

SLIDER
17
469
269
502
testing-rate-symptomatic
testing-rate-symptomatic
0
1
0.1
0.001
1
NIL
HORIZONTAL

SLIDER
14
319
154
352
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
163
317
310
350
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
14
277
151
310
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
10
369
156
402
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
161
367
311
400
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
770
587
943
620
initial-exposed
initial-exposed
0
1000
500.0
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
3.0
1
1
NIL
HORIZONTAL

SLIDER
333
326
506
359
flow-rate
flow-rate
0
1
0.0
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
334
29
507
62
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
338
10
422
28
Population\n
12
0.0
1

TEXTBOX
19
6
103
24
Pandemic
12
0.0
1

TEXTBOX
14
258
98
276
Mortality
12
0.0
1

TEXTBOX
334
306
418
324
Connectivity
12
0.0
1

TEXTBOX
17
449
148
479
Control and testing
12
0.0
1

SLIDER
334
104
506
137
pop-sd-multiplier
pop-sd-multiplier
0.01
1.2
0.62
0.01
1
NIL
HORIZONTAL

MONITOR
336
253
445
298
total-pop
sum [pop-0] of locales
0
1
11

MONITOR
334
200
442
245
max-pop
max [pop-0] of locales
0
1
11

MONITOR
333
150
413
195
min-pop
min [pop-0] of locales
0
1
11

INPUTBOX
349
457
498
517
alert-levels-R0
[2.5 2.1 1.7 1.3 0.9]
1
0
String

MONITOR
20
26
92
72
mean-R0
mean-R0
3
1
11

SLIDER
333
540
506
574
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
1460
19
1518
65
dead
sum [dead] of locales
0
1
11

SWITCH
770
628
917
662
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
18
508
295
542
testing-rate-presymptomatic
testing-rate-presymptomatic
0
0.1
0.005
0.001
1
NIL
HORIZONTAL

SLIDER
17
547
237
581
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
log (tests + 1)
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"tests" 1.0 0 -16777216 true "" "plot log (total-tests + 1) 10"
"tests-pos" 1.0 0 -7500403 true "" "plot log (total-tests-positive + 1) 10"
"new-tests" 1.0 0 -2674135 true "" "plot log (total-new-tests + 1) 10"
"new-tests-pos" 1.0 0 -955883 true "" "plot log (total-new-tests-positive + 1) 10"

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
NetLogo 6.1.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
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
