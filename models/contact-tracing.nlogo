breed [cases case]

cases-own [
  time-0
  base-R
  current-R
  clinical?
  isolated?
  recovered?

  exposures ;; continously maintained list of exposures from this case

  t-onset-1
  control

  n-initial
]


to setup
  clear-all
  reset-ticks
  ask patches [ set pcolor white ]

  if use-seed? [ random-seed seed ]

  set-default-shape cases "circle"

  create-cases 1 [
    initialise-case ticks p-clinical nobody
  ]
end


to initialise-case [t p-clin parent]
  if debug? [ show word "New case at t " t ]
  if parent != nobody [
    create-link-from parent
  ]
  set time-0 t
  set clinical? (random-float 1 < p-clin)        ;; TPM paper
  ifelse clinical? [
    set base-R R-clinical
    set color [255 0 0 255]
  ]
  [
    set base-R R-clinical * subclin-fraction-R
    set color [0 0 255 255]
  ]
  set exposures get-exposures
  let all-new-exposures split-list-at-value exposures t
  set exposures last all-new-exposures

  let newer-exposures first all-new-exposures
  add-cases-from-list newer-exposures

  set t-onset-1 random-gamma 5.5 0.95             ;; TPM paper

  set isolated? false
  set recovered? false
end


to go
  add-cases (ticks + 1)
  trace-cases (ticks + 1)
  progress-cases
  layout-radial cases links case 0
  draw-cases
  tick
end


to add-cases [time-now]
  ask cases with [not recovered?] [
    let split-exposures split-list-at-value exposures time-now
    let new-cases first split-exposures
    set exposures last split-exposures

    if debug? [show word "Adding cases " new-cases]
    add-cases-from-list new-cases
  ]
end


to add-cases-from-list [L]
  foreach L [ t ->
    add-case t self
  ]
end


to add-case [t parent-case]
  ask parent-case [
    if random-float 1 < control [
      hatch 1 [
        initialise-case t p-clinical myself
      ]
    ]
  ]
end


to-report get-exposures
  set n-initial random-poisson get-lambda
  report sort n-values n-initial [t ->
    time-0 + random-weibull w-shape weibull-scale-given-shape
  ]
end

to-report get-lambda
  report base-R
end


to draw-cases
  ask cases [
    if isolated? [
      set color grey
    ]
    if recovered? [
      set shape "x"
      set color green
    ]
  ]
end


to progress-cases
  ask cases with [not recovered?] [
    set recovered? length exposures = 0
  ]
  ask cases [
    set control lockdown-control * ifelse-value isolated? [isolation-control] [1]
  ]
end


to trace-cases [time-now]
  ;; since only clinical cases attract attention...
  ask cases with [clinical? and not isolated? and time-now > (time-0 + t-onset-1)] [
    let time-since-onset time-now - time-0 - t-onset-1
    set isolated? random-float 1 < p-detection time-since-onset
  ]
  ask cases with [isolated?] [
    predecessor-trace
    successor-trace
  ]
end

to predecessor-trace
  ask in-link-neighbors [
    set isolated? clinical? and random-float 1 < p-contact-detection
    if isolated? [
      predecessor-trace
      successor-trace
    ]
  ]
end

to successor-trace
  ask out-link-neighbors [
    set isolated? clinical? and random-float 1 < p-contact-detection
    if isolated? [
      successor-trace
    ]
  ]
end

to-report p-detection [t]
  report 1 - (1 - p-detection-per-day) ^ t
end


;; ---------------------------------------
;; SOME STRING and LIST HANDLING STUFF
;; ---------------------------------------
to-report insert-value-in-order [L x]
  report (sentence filter [v -> v < x] L x filter [v -> v >= x] L)
end

to-report insert-values-in-order [L new]
  report reduce [ [a b] -> insert-value-in-order a b] (fput L new)
end

to-report split-list-at-value [L x] ;; assume list is ordered
  report (list filter [v -> v < x] L filter [v -> v >= x] L)
end



;; -----------------------------------
;; HANDY DISTRIBUTIONS!
;; -----------------------------------

; This binomial algorithm from
; Devroye. L. 1960. Generating the maximum of independent identically
; distributed random variables. Computers and Mathematics with
; Applications 6, 305-315.
; should be a bit quicker because it only needs ~ np random-float calls
to-report random-binomial [n p]
  if p = 0 [ report 0 ]
  if p = 1 [ report n ]
  let y 0
  let x 0
  ; since O(np) use fact that Bin(n, p) = n - Bin(n, 1-p)
  ifelse p < 0.5 [
    let c ln (1 - p)
    loop [
      set y y + int (ln (random-float 1) / c) + 1
      ifelse y < n
      [ set x x + 1 ]
      [ report x ]
    ]
  ]
  [
    let c ln p
    loop [
      set y y + int (ln (random-float 1) / c) + 1
      ifelse y < n
      [ set x x + 1 ]
      [ report n - x ]
    ]
  ]
end

;; https://stackoverflow.com/questions/33611708/random-number-generator-with-generalized-pareto-distribution-and-weilbull-distri
to-report random-weibull [shp scale]
  report scale * (-1 * ln (random-float 1)) ^ (1 / shp)
end

to-report cumulative-weibull [shp scale x]
  if x <= 0 [
    report 0
  ]
  report 1 - 1 / exp ((x / scale) ^ shp)
end

to-report weibull-in-interval [shp scale t1-t2]
  report cumulative-weibull shp scale (item 1 t1-t2) - cumulative-weibull shp scale (item 0 t1-t2)
end

to-report dweibull [x shp scale]
  if x <= 0 [ report 0 ]
  report (shp / scale) * ((x / scale) ^ (shp - 1)) / exp ((x / scale) ^ shp)
end

to-report weibull-scale-given-shape
  report mode * (w-shape / (w-shape - 1)) ^ (1 / w-shape)
end
@#$#@#$#@
GRAPHICS-WINDOW
210
10
626
427
-1
-1
8.0
1
10
1
1
1
0
0
0
1
-25
25
-25
25
1
1
1
days
30.0

BUTTON
57
12
122
47
NIL
setup\n\n
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
130
14
194
48
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
14
190
199
223
R-clinical
R-clinical
0
10
3.0
0.001
1
NIL
HORIZONTAL

SLIDER
13
228
200
261
subclin-fraction-R
subclin-fraction-R
0
1
0.5
0.001
1
NIL
HORIZONTAL

SLIDER
12
264
202
298
p-clinical
p-clinical
0
1
0.667
0.001
1
NIL
HORIZONTAL

SLIDER
12
365
184
398
mode
mode
0
10
4.86
0.01
1
NIL
HORIZONTAL

SLIDER
12
402
184
435
w-shape
w-shape
0
5
2.833
0.001
1
NIL
HORIZONTAL

MONITOR
15
442
78
487
w-scale
weibull-scale-given-shape
4
1
11

BUTTON
132
60
197
94
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

SLIDER
640
273
824
306
lockdown-control
lockdown-control
0
1
0.72
0.01
1
NIL
HORIZONTAL

SLIDER
20
104
193
137
seed
seed
0
100
51.0
1
1
NIL
HORIZONTAL

BUTTON
39
60
122
94
go-10
repeat 10 [go]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SWITCH
84
141
194
174
use-seed?
use-seed?
1
1
-1000

SWITCH
768
492
886
525
debug?
debug?
1
1
-1000

SLIDER
641
57
854
90
p-detection-per-day
p-detection-per-day
0
1
0.17
0.01
1
NIL
HORIZONTAL

SLIDER
640
162
853
195
p-contact-detection
p-contact-detection
0
1
0.9
0.01
1
NIL
HORIZONTAL

SLIDER
643
381
827
415
isolation-control
isolation-control
0
1
0.65
0.01
1
NIL
HORIZONTAL

MONITOR
219
440
313
485
NIL
count cases
0
1
11

TEXTBOX
14
330
199
360
Mode is peak infectiousness \npost exposure
12
0.0
1

TEXTBOX
646
24
863
54
Probability of detection per day \nafter onset of symptoms
12
0.0
1

TEXTBOX
642
109
880
156
Probability of detection of contacts \nboth back and forward in the chain \n(stopped by subclinical cases)
12
0.0
1

TEXTBOX
645
210
837
267
general lowering of R due to \nsocial distancing measures: \n0.32 for L4, 0.52 for L3, \n0.72 for L2 and 1 for L1
12
0.0
1

TEXTBOX
644
332
839
379
Additional damping of R due\nto isolation of contact traced\ncases
12
0.0
1

TEXTBOX
326
437
628
527
Green crosses are recovered cases, which is occurs when they've finished exposing cases in their queue. Greys have been isolated by contact tracing. Reds are clinical, blues subclinical.
12
0.0
1

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
