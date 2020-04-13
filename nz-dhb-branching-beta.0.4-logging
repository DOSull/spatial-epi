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
  cum-cases          ;; exposed
  cum-infected
  cum-recovered      ;; recovered

  new-cases
  new-infected
  new-recovered

  recent-new-cases   ;; list of recent cases per 'time-horizon' setting
  recent-tests

  exposures          ;; list of times of exposures 'in the queue'
  new-infections     ;; new infections this tick

  alert-level        ;; local alert level which controls...
  my-alert-indicator ;; local level turtle indicator
  my-control         ;; control level

  name
]

breed [cases case]
cases-own [
  base-R
  time-0
  t-onset-1
  t-onset-2
  clinical?
  isolated?
  in-hospital?
  time-leaving-hospital
  my-locale
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
  control-levels
  flow-levels
  trigger-levels

  total-cases          ;; cases
  total-infected       ;; cases
  total-recovered      ;; recovered

  log-header-file-name
  log-file-name
  full-file-name
  date-time
  model-name
  labels-on?

  alert-level-changes
]

to setup
  clear-all
  reset-ticks

  ask patches [set pcolor cyan + 1]

  set-default-shape locales "circle"
  set-default-shape cases "circle"
  set-default-shape levels "square 3"
  set-default-shape connections "myshape"

  if use-seed? [random-seed seed]

  set control-levels initialise-control-levels
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

  let size-adjust 0.25 / sqrt (count locales / count patches)
  ask (turtle-set locales levels) [ set size size * size-adjust ]
  ask connections [ set thickness thickness * size-adjust ]

  ;; initial exposures
  ifelse uniform-by-pop? [
    uniform-cases-by-population initial-infected
  ]
  [
    repeat initial-infected [
      ask one-of locales [
        initial-exposure (0 - random-exponential 5)
      ]
    ]
  ]
  let pre-start-exposed-locales locales with [length exposures > 0]
  if any? pre-start-exposed-locales [
    let pre-t floor min [item 0 exposures] of pre-start-exposed-locales
    while [pre-t < 0] [
      set pre-t pre-t + 1
      run-one-day pre-t
    ]
  ]
  update-plots

  if log-all-locales? [ initialise-logging ]

  set labels-on? true
  redraw
  paint-land lime - 1 4
end

;; add an exposure to the queue
to initial-exposure [t]
  set exposures insert-value-in-order exposures (t + random-weibull 2.83 5.67)
end


;; susceptible weighted infection so
;; that higher population locales are more
;; likely to receive exposures
to uniform-cases-by-population [n]
  repeat n [
    ask random-locale-by-pop [
      initial-exposure (0 - random-exponential 5)
    ]
  ]
end

to-report random-locale-by-pop
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
  report item idx ordered-locales
end

to-report cumulative-sum [lst]
  report but-first reduce [ [a b] -> lput (last a + b) a] (fput [0] lst)
end


to update-parameters
  set total-cases count cases with [clinical?]
  set total-infected count cases
  set total-recovered total-recovered + sum [new-recovered] of locales

  ask locales [
    set cum-cases cum-cases + new-cases
    set cum-infected cum-infected + new-infected
    set cum-recovered cum-recovered + new-recovered

    set recent-new-cases fput new-cases recent-new-cases

    let new-tests new-cases + (random-binomial (new-infected - new-cases) test-rate-presymp) + (random-binomial susceptible test-rate-gen)
    set recent-tests fput new-tests recent-tests
  ]
end


;; -------------------------
;; Main
;; -------------------------
to go
  if (sum [length exposures] of locales + count cases = 0 and ticks > 0) [
    update-parameters
    update-plots
    stop
  ]

  repeat random-poisson new-exposures-arriving [
    ask random-locale-by-pop [
      initial-exposure ticks
    ]
  ]

  run-one-day (ticks + 1)

  if log-all-locales? [
    file-open full-file-name
    ask locales [
      file-print output-locale
    ]
    file-close
  ]

  redraw

  if ticks >= start-lifting-quarantine and (ticks - start-lifting-quarantine) mod time-horizon = 0 [
    change-alert-levels
  ]

  tick
end


to run-one-day [time]
  flush-exposures time
  ask cases [
    infect-others time
  ]
  ask cases [
    progress-infection time
  ]

  update-parameters
end


to flush-exposures [time-now]
  ask locales [
    let split-exposures split-list-at-value exposures time-now
    set new-infections item 0 split-exposures
    set exposures item 1 split-exposures

    set new-cases 0
    set new-infected 0
    set new-recovered 0
    while [ length new-infections > 0 ] [
      let t first new-infections
      ;; remove it immediately from queue
      ;; in case a new one sneaks in in front of it!
      set new-infections but-first new-infections
      hatch 1 [
        initialise-case t
      ]
      set new-infected new-infected + 1
      set susceptible susceptible - 1
    ]
  ]
end

;; ------------------------------
;; ALERT LEVELS
;; ------------------------------
to change-alert-levels
  if alert-policy = "static" [
    ask locales [
      let new-level initial-alert-level
      if new-level != alert-level [
        set alert-level new-level
        set alert-level-changes alert-level-changes + 1
      ]
    ]
    enact-new-levels
    stop
  ]

  if alert-policy = "local-random" [
    ask locales [
      let alert-level-change one-of (list -1 0 1)
      let new-level clamp (alert-level + alert-level-change) 1 4
      if new-level != alert-level [
        set alert-level new-level
        set alert-level-changes alert-level-changes + 1
      ]
    ]
    enact-new-levels
    stop
  ]

  if alert-policy = "local" [
    ask locales [
      let current-new-cases recent-total recent-new-cases
      let current-tests recent-total recent-tests
      if current-tests > 0 [
        let positive-test-rate current-new-cases / current-tests
        let new-level get-new-alert-level positive-test-rate
        if new-level != alert-level [
          ifelse new-level < alert-level [
            set alert-level clamp (alert-level - 1) 1 4
          ]
          [
            set alert-level new-level
          ]
          set alert-level-changes alert-level-changes + 1
        ]
      ]
    ]
    enact-new-levels
    stop
  ]

  if alert-policy = "global" [
    let current-new-cases sum [recent-total recent-new-cases] of locales
    let current-tests sum [recent-total recent-tests] of locales
    if current-tests > 0 [
      let positive-test-rate current-new-cases / current-tests
      let new-level get-new-alert-level positive-test-rate
      let a [alert-level] of one-of locales
      if new-level != a [
        ifelse new-level < a [
          ask locales [
            set alert-level clamp (alert-level - 1) 1 4
          ]
        ]
        [
          ask locales [
            set alert-level new-level
          ]
        ]
        set alert-level-changes alert-level-changes + count locales
      ]
    ]
    enact-new-levels
    stop
  ]
end

to-report recent-total [lst]
  report sum sublist lst 0 time-horizon
end

to-report get-new-alert-level [rate]
  report 1 + length filter [x -> rate > x] trigger-levels
end

to-report clamp [x mn mx]
  report max (list mn min (list mx x))
end

to enact-new-levels
  ask locales [
    set my-control item (alert-level - 1) control-levels
  ]
  ask levels [
    set alert-level [alert-level - 1] of my-locale
    draw-level
  ]
  ask connections [
    set my-flow-rate min [item (alert-level - 1) flow-levels] of both-ends
  ]
end



; ----------------------------------------------
; BRANCHING MODEL GROWTH
; ----------------------------------------------
to infect-others [time-now]
  let time-since-infection time-now - time-0

  let R base-R
  if clinical? and not isolated? and time-since-infection > (t-onset-1 + t-onset-2) [
    set isolated? true
    set base-R base-R * 0.65
  ]

  let n random-poisson get-lambda time-since-infection
  let new-infection-times split-list-at-value sort n-values n [t -> time-0 + random-weibull 2.83 5.67] time-now

  let the-locale my-locale
  if not clinical?  [
    ask my-locale [
      if random-float 1 < (item (alert-level - 1) flow-levels) [
        set the-locale choose-neighbour-by-weight
      ]
    ]
  ]
  ask the-locale [
    set new-infections insert-values-in-order new-infections (item 0 new-infection-times)
    set exposures insert-values-in-order exposures (item 1 new-infection-times)
  ]
end


to-report get-lambda [t]
  let R base-R * [my-control] of my-locale
  let lambda R * [susceptible / pop-0] of my-locale * weibull-in-interval 2.83 5.67 (t - 1) t
  if t > (t-onset-1 + t-onset-2) [
    set lambda lambda * 0.65 * [my-control] of my-locale
  ]
  report lambda
end

to make-new-cases [n t]
  ask my-locale [
    hatch n [ initialise-case t ]
  ]
end


to progress-infection [time-now]
  let time-since-infection time-now - time-0

  if time-since-infection > 30 [
;    ifelse in-hospital? [
;      ask my-locale [
;        set new-dead new-dead + 1
;        set pop-0 pop-0 - 1
;      ]
;    ]
;    [
      ask my-locale [
        set new-recovered new-recovered + 1
;     ]
    ]
    die
  ]
  if time-now = ceiling (time-0 + t-onset-1) [ ;; I think we only check for hospitalisation once
    set in-hospital? random-float 1 < 0.078
    set time-leaving-hospital t-onset-1 + random-exponential 10
  ]
  if in-hospital? [
    set in-hospital? time-now > time-leaving-hospital
  ]
end


;to-report get-effective-infected
;  report infected + flow-rate * sum [my-flow-rate * w * [infected] of other-end] of my-in-connections
;end



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
    let xy random-xy-with-buffer 0.1
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
  report (list rescale random-xcor (min-pxcor - 0.5) (max-pxcor + 0.5) buffer true
               rescale random-ycor (min-pycor - 0.5) (max-pycor + 0.5) buffer false)
end

to initialise-locales
  let mean-pop-0 mean [pop-0] of locales
  ask locales [
    set size (pop-0 / mean-pop-0) ^ (1 / 3)
    set susceptible pop-0
    set new-cases 0
    set new-recovered 0

    set cum-cases 0
    set cum-recovered 0

    set recent-new-cases []
    set recent-tests []

    set exposures []

    set alert-level initial-alert-level
    set my-control item (alert-level - 1) control-levels
  ]
end


to initialise-locales-from-string [s]
  let locales-data but-first split-string s "\n"

  let xs map [x -> read-from-string item 2 split-string x " "] locales-data
  let ys map [y -> read-from-string item 3 split-string y " "] locales-data
  let min-x min xs    let max-x max xs
  let min-y min ys    let max-y max ys
  set xs map [x -> rescale x min-x max-x 0.1 true] xs
  set ys map [y -> rescale y min-y max-y 0.1 false] ys

  create-locales length locales-data

  (foreach locales-data sort locales xs ys [ [line loc x y] ->
    let parameters split-string line " "
    ask loc [
      set name item 1 parameters
      set label name
      set label-color black
      set pop-0 read-from-string item 4 parameters
      setxy x y
    ]
  ])
  initialise-locales
end

to-report rescale [z min-z max-z buffer x?]
  let new-range ifelse-value x? [(world-width - 1) * (1 - buffer)] [(world-height - 1) * (1 - buffer)]
  let new-min ifelse-value x? [(world-width - 1 - new-range) / 2] [(world-height - 1 - new-range) / 2]
  let new-z new-min + (z - min-z) / (max-z - min-z) * new-range
  report new-z
end


to setup-levels
  ask locales [
    let x nobody
    hatch 1 [
      set size size * 1.2
      set breed levels
      set my-locale myself
      set x self
      set label ""
    ]
    set my-alert-indicator x
  ]
end


to initialise-case [t]
  set breed cases
  set my-locale myself

  set time-0 t
  set clinical? (random-float 1 < p-clinical)        ;; TPM paper
  ifelse clinical? [
    set base-R R-clin
    set color red
    ask my-locale [
      set new-cases new-cases + 1
    ]
  ]
  [
    set base-R R-clin * 0.5
    set color blue
  ]

  set t-onset-1 random-gamma 5.5 0.95             ;; TPM paper
  ifelse fast-isolation? [
    set t-onset-2 random-exponential 2.18     ;; TPM paper
  ]
  [
    set t-onset-2 random-exponential 6
  ]

  set isolated? false
  set in-hospital? false

  set size 0.1
  set label ""
  set heading random 360
  jump random-float [size / 2] of my-locale
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


to initialise-connections-from-string [s]
  let edges but-first split-string s "\n"
  foreach edges [ edge ->
    let parameters split-string edge " "
    let v1 locale (read-from-string item 0 parameters)
    let v2 locale (read-from-string item 1 parameters)
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
  set color [0 102 204 127]
  set w 1 / d
  set thickness w
end

to reweight-connections
  ask locales [
    let total-w sum [w] of my-out-connections
    let w-correction 1 / total-w
    let thickness-correction 1 / total-w
    ask my-out-connections [
      set w w * w-correction
      set thickness thickness * thickness-correction
    ]
  ]
end

to-report choose-neighbour-by-weight
  let neighbours reverse sort-on [w] my-out-connections
  let cumulative-weights cumulative-sum map [x -> [w] of x] neighbours
  let picker random-float last cumulative-weights
  ;show map [x -> precision x 4] cumulative-weights
  ;show picker
  let idx length filter [x -> x <= picker] cumulative-weights
  report [other-end] of item idx neighbours
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
  set color scale-color red (count cases with [my-locale = myself] / pop-0) pop-0 0
end

to toggle-labels
  set labels-on? not labels-on?
  ask locales [
    ifelse labels-on?
    [ set label name ]
    [ set label "" ]
  ]
end

to draw-level
  set alert-level [alert-level] of my-locale
  set color item alert-level (list pcolor lime yellow orange red)
end


;; --------------------------------------------------
;; COLOUR IN SOME 'LAND'
;; --------------------------------------------------
to paint-land [c n]
  let loc nobody
  let t nobody
  ask locales [
    set loc self
    ask patch-here [
      sprout 1 [
        set color c
        move-to loc
        set pcolor color
        set t self
      ]
    ]
    let targets sort-on [length-link] my-out-connections
    set targets sublist targets 0 min (list length targets n)
    foreach targets [ edge ->
      ask t [
        walk-edge self loc edge
      ]
    ]
    ask t [die]
  ]
end

to walk-edge [ttl loc edge]
  ask loc [
    let tgt [other-end] of edge
    ask ttl [
      face tgt
      let d distance tgt
      let ceiling-d ceiling d
      repeat ceiling-d [
        fd d / ceiling-d
        set pcolor color
      ]
    ]
  ]
end

to-report length-link
  let d 0
  ask one-of both-ends [
    set d distance [other-end] of myself
  ]
  report d
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

to-report weibull-in-interval [shp scale x1 x2]
  report cumulative-weibull shp scale x2 - cumulative-weibull shp scale x1
end

;; ---------------------------------------
;; SOME STRING and LIST HANDLING STUFF
;; ---------------------------------------
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


to-report insert-value-in-order [L x]
  report (sentence filter [v -> v < x] L x filter [v -> v >= x] L)
end

to-report insert-values-in-order [L new]
  report reduce [ [a b] -> insert-value-in-order a b] (fput L new)
end

to-report split-list-at-value [L x] ;; assume list is ordered
  report (list filter [v -> v < x] L filter [v -> v >= x] L)
end


to-report replace [s a b]
  let i position a s
  report replace-item i s b
end

to-report join-list [lst sep]
  report reduce [ [a b] -> (word a sep b) ] lst
end

to-report initialise-control-levels
  let control-level-specs (split-string alert-levels-control "\n")
  let control-level-names map [spec -> item 0 split-string spec " "] control-level-specs
  let idx position control-scenario control-level-names
  report read-from-string (join-list (but-first split-string item idx control-level-specs " ") " ")
end

;; ----------------------------------------
;; LOGGING
;; ----------------------------------------
to initialise-logging
  set model-name "nz-dhb-branching-beta.0.3-logging"
  set date-time replace date-and-time "." ":"
  let base-file-name (word model-name "-" date-time)
  set log-file-name (word base-file-name ".csv")
  set full-file-name (word log-folder "/" log-file-name)
  set log-header-file-name (word log-folder "/" base-file-name ".header")
  if log-all-locales? [
    file-open log-header-file-name
    file-print log-file-header
    file-close
    file-open full-file-name
    file-print output-locale-header
    file-close
  ]
end

to-report output-locale
  report join-list (list ticks who name pop-0 susceptible cum-cases cum-infected cum-recovered new-cases new-infected new-recovered alert-level) ","
end

to-report output-locale-header
  report "ticks,who,name,pop.0,susceptible,cum.cases,cum.infected,cum.recovered,new.cases,new.infected,new.recovered,alert.level"
end

to-report log-file-header
  let parameters (list "name,value")
  set parameters lput join-list (list "model.name" model-name) "," parameters
  set parameters lput join-list (list "log.file.name" log-file-name) "," parameters
  set parameters lput join-list (list "log.all.locales?" log-all-locales?) "," parameters
  set parameters lput join-list (list "log.folder" log-folder) "," parameters

  set parameters lput join-list (list "use.seed?" use-seed?) "," parameters
  set parameters lput join-list (list "seed" seed) "," parameters

  set parameters lput join-list (list "initial.infected" initial-infected) "," parameters
  set parameters lput join-list (list "uniform.by.pop?" uniform-by-pop?) "," parameters
  set parameters lput join-list (list "initial.alert.level" initial-alert-level) "," parameters
  set parameters lput join-list (list "alert.policy" alert-policy) "," parameters
  set parameters lput join-list (list "start.lifting.quarantine" start-lifting-quarantine) "," parameters
  set parameters lput join-list (list "time.horizon" time-horizon) "," parameters

  set parameters lput join-list (list "alert.levels.flow" alert-levels-flow) "," parameters

  set parameters lput join-list (list "test.rate.symp" test-rate-symp) "," parameters
  set parameters lput join-list (list "test.rate.presymp" test-rate-presymp) "," parameters
  set parameters lput join-list (list "test.rate.gen" test-rate-gen) "," parameters

  set parameters lput join-list (list "population" population) "," parameters
  set parameters lput join-list (list "num.locales" num-locales) "," parameters
  set parameters lput join-list (list "pop.sd.multiplier" pop-sd-multiplier) "," parameters

  set parameters lput join-list (list "trigger.levels" trigger-levels) "," parameters
  set parameters lput join-list (list "flow.levels" flow-levels) "," parameters
  set parameters lput join-list (list "control.levels" control-levels) "," parameters
  set parameters lput join-list (list "control.scenario" control-scenario) "," parameters
  set parameters lput join-list (list "fast.isolation" fast-isolation?) "," parameters

  set parameters lput join-list (list "R.clin" r-clin) "," parameters
  set parameters lput join-list (list "p.clinical" p-clinical) "," parameters

  set parameters lput join-list (list "min.pxcor" min-pxcor) "," parameters
  set parameters lput join-list (list "max.pxcor" max-pxcor) "," parameters
  set parameters lput join-list (list "min.pycor" min-pycor) "," parameters
  set parameters lput join-list (list "max.pycor" max-pycor) "," parameters

  report join-list parameters "\n"
end
@#$#@#$#@
GRAPHICS-WINDOW
540
12
1028
741
-1
-1
30.0
1
14
1
1
1
0
0
0
1
0
15
0
23
1
1
1
days
100.0

BUTTON
457
17
531
51
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
458
175
532
210
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
458
62
532
97
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
9
62
183
95
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
10
693
170
726
test-rate-symp
test-rate-symp
0
1
0.5
0.01
1
NIL
HORIZONTAL

SLIDER
380
225
529
258
initial-infected
initial-infected
0
2000
1310.0
10
1
NIL
HORIZONTAL

SLIDER
264
58
437
91
seed
seed
0
100
30.0
1
1
NIL
HORIZONTAL

SWITCH
298
19
437
52
use-seed?
use-seed?
1
1
-1000

SLIDER
9
25
182
58
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
14
7
98
25
Population\n
12
0.0
1

TEXTBOX
1039
36
1123
54
Pandemic
12
0.0
1

TEXTBOX
10
314
94
332
Connectivity
12
0.0
1

TEXTBOX
12
451
144
470
Control and testing
12
0.0
1

SLIDER
9
100
181
133
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
99
144
179
189
max-pop
max [pop-0] of locales
0
1
11

MONITOR
11
145
92
190
min-pop
min [pop-0] of locales
0
1
11

SLIDER
375
343
529
376
initial-alert-level
initial-alert-level
1
4
4.0
1
1
NIL
HORIZONTAL

SWITCH
380
267
527
300
uniform-by-pop?
uniform-by-pop?
0
1
-1000

PLOT
1042
154
1402
437
Cumulative totals
Days
Log (count + 1}
0.0
10.0
0.0
3.0
true
true
"" ""
PENS
"cases" 1.0 0 -2674135 true "" "plot log (sum [cum-cases] of locales + 1) 10"
"infected" 1.0 0 -13345367 true "" "plot log (sum [cum-infected] of locales + 1) 10"
"recovered" 1.0 0 -13840069 true "" "plot log (sum [cum-recovered] of locales + 1) 10"

SLIDER
10
733
170
766
test-rate-presymp
test-rate-presymp
0
0.1
0.05
0.001
1
NIL
HORIZONTAL

SLIDER
9
773
169
806
test-rate-gen
test-rate-gen
0
0.002
0.001
0.0001
1
NIL
HORIZONTAL

INPUTBOX
313
435
532
495
alert-levels-flow
[0.1 0.05 0.025 0.01]
1
0
String

CHOOSER
375
382
528
427
alert-policy
alert-policy
"global" "local" "local-random" "static"
0

MONITOR
1045
100
1119
145
infected
count cases
0
1
11

MONITOR
1320
101
1393
146
recovered
total-recovered
0
1
11

INPUTBOX
313
499
531
559
alert-level-triggers
[0.0001 0.00025 0.0005 1]
1
0
String

SLIDER
329
603
533
636
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
548
752
719
785
log-all-locales?
log-all-locales?
1
1
-1000

SLIDER
329
565
531
598
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
457
320
536
338
Alert levels
12
0.0
1

MONITOR
220
647
278
692
pop-lev-0
sum [pop-0] of locales with [alert-level = 0]
0
1
11

MONITOR
284
647
342
692
pop-lev-1
sum [pop-0] of locales with [alert-level = 1]
0
1
11

MONITOR
346
647
404
692
pop-lev-2
sum [pop-0] of locales with [alert-level = 2]
0
1
11

MONITOR
409
647
468
692
pop-lev-3
sum [pop-0] of locales with [alert-level = 3]
0
1
11

MONITOR
474
647
533
692
pop-lev-4
sum [pop-0] of locales with [alert-level = 4]
0
1
11

MONITOR
220
698
278
743
n-lev-0
count locales with [alert-level = 0]
0
1
11

MONITOR
284
698
342
743
n-lev-1
count locales with [alert-level = 1]
0
1
11

MONITOR
346
698
404
743
n-lev-2
count locales with [alert-level = 2]
0
1
11

MONITOR
409
698
468
743
n-lev-3
count locales with [alert-level = 3]
0
1
11

MONITOR
474
698
533
743
n-lev-4
count locales with [alert-level = 4]
0
1
11

INPUTBOX
1043
746
1228
806
log-folder
test
1
0
String

SWITCH
204
165
390
198
initialise-from-nz-data?
initialise-from-nz-data?
0
1
-1000

INPUTBOX
10
210
339
305
dhbs
ID name x y pop\n0 Northland 1640809 6107096 178500\n1 Waitemata 1704007 6060146 622350\n2 Auckland 1794796 5977673 538840\n3 Counties.Manukau 1752361 5844360 557790\n4 Waikato 1868112 5737176 417130\n5 Lakes 1897334 5643762 110040\n6 Bay.of.Plenty 1981173 5745420 236820\n7 Tairawhiti 2043079 5706142 48965\n8 Taranaki 1713209 5586361 119640\n9 Hawkes.Bay 1954114 5586840 165360\n10 Whanganui 1821171 5588135 64565\n11 MidCentral 1825104 5499583 178240\n12 Hutt.Valley 1790801 5403979 149270\n13 Capital.and.Coast 1727799 5412778 316620\n14 Wairarapa 1854566 5434609 44845\n15 Nelson.Marlborough 1581787 5353746 150280\n16 West.Coast 1466785 5271651 32445\n17 Canterbury 1493549 5180880 562870\n18 South.Canterbury 1352314 5083543 60025\n19 Southern 1279834 4972310 328230\n
1
1
String

INPUTBOX
10
332
338
425
connectivity
ID1 ID2 cost weight weight2\n8 9 0.227069 4.40394 19.394684\n8 13 0.193776 5.160594 26.631733\n8 10 0.093423 10.704056 114.576815\n8 7 0.345958 2.890528 8.355151\n8 3 0.198042 5.04943 25.496742\n8 0 0.335123 2.983978 8.904126\n8 5 0.167174 5.98179 35.78181\n8 6 0.199933 5.001679 25.016795\n8 1 0.236637 4.225886 17.858116\n8 2 0.217664 4.594228 21.106932\n8 4 0.147367 6.785774 46.046732\n8 11 0.131898 7.581596 57.480603\n8 12 0.196719 5.083381 25.840766\n8 14 0.1913 5.227401 27.325723\n8 19 0.702723 1.423036 2.025031\n8 18 0.546829 1.828725 3.344233\n8 17 0.452491 2.209987 4.884042\n8 15 0.309359 3.23249 10.448993\n8 16 0.465441 2.148499 4.616047\n9 8 0.227069 4.40394 19.394684\n9 13 0.174988 5.714684 32.657614\n9 10 0.143124 6.986926 48.817133\n9 7 0.154221 6.484195 42.044779\n9 3 0.229013 4.366559 19.066835\n9 0 0.366094 2.731537 7.461295\n9 5 0.119525 8.366462 69.997687\n9 6 0.164249 6.088325 37.067702\n9 1 0.267608 3.736811 13.963754\n9 2 0.248636 4.021951 16.176089\n9 4 0.171767 5.821835 33.893763\n9 11 0.103927 9.622145 92.585678\n9 12 0.17682 5.655482 31.984481\n9 14 0.134911 7.412309 54.94233\n9 19 0.683935 1.462128 2.137818\n9 18 0.528041 1.893793 3.586451\n9 17 0.433703 2.305725 5.316368\n9 15 0.290571 3.441503 11.843945\n9 16 0.446653 2.238875 5.01256\n13 8 0.193776 5.160594 26.631733\n13 9 0.174988 5.714684 32.657614\n13 10 0.109831 9.104883 82.898897\n13 7 0.326629 3.061574 9.373238\n13 3 0.353972 2.825081 7.98108\n13 0 0.491053 2.036439 4.147086\n13 5 0.246502 4.056755 16.457265\n13 6 0.291226 3.433755 11.790674\n13 1 0.392567 2.547337 6.488927\n13 2 0.373594 2.676699 7.164719\n13 4 0.297856 3.357331 11.271675\n13 11 0.07454 13.415666 179.980082\n13 12 0.012081 82.771665 6851.148445\n13 14 0.053941 18.538685 343.68286\n13 19 0.509581 1.962397 3.851002\n13 18 0.353687 2.827358 7.993953\n13 17 0.259349 3.855804 14.867228\n13 15 0.116217 8.604601 74.039161\n13 16 0.272299 3.672432 13.486755\n10 8 0.093423 10.704056 114.576815\n10 9 0.143124 6.986926 48.817133\n10 13 0.109831 9.104883 82.898897\n10 7 0.294766 3.392522 11.509203\n10 3 0.257996 3.876026 15.023576\n10 0 0.395077 2.531151 6.406726\n10 5 0.164657 6.073228 36.884097\n10 6 0.209381 4.775982 22.810001\n10 1 0.296591 3.371649 11.368015\n10 2 0.277618 3.602066 12.974877\n10 4 0.207321 4.823434 23.265513\n10 11 0.047953 20.853585 434.871996\n10 12 0.112774 8.867253 78.628178\n10 14 0.107355 9.314917 86.767676\n10 19 0.618778 1.616088 2.611742\n10 18 0.462884 2.160367 4.667187\n10 17 0.368546 2.713362 7.362335\n10 15 0.225414 4.43628 19.680584\n10 16 0.381496 2.621257 6.87099\n7 8 0.345958 2.890528 8.355151\n7 9 0.154221 6.484195 42.044779\n7 13 0.326629 3.061574 9.373238\n7 10 0.294766 3.392522 11.509203\n7 3 0.270829 3.692362 13.633541\n7 0 0.40778 2.452302 6.013784\n7 5 0.182956 5.465785 29.874806\n7 6 0.157261 6.358875 40.435291\n7 1 0.309294 3.233172 10.453399\n7 2 0.290321 3.444457 11.864286\n7 4 0.225684 4.430975 19.633541\n7 11 0.255568 3.912846 15.310363\n7 12 0.328461 3.0445 9.268983\n7 14 0.286552 3.489765 12.178459\n7 19 0.835576 1.196779 1.43228\n7 18 0.679682 1.471275 2.164651\n7 17 0.585345 1.708395 2.918615\n7 15 0.442212 2.261358 5.113739\n7 16 0.598294 1.671418 2.793637\n3 8 0.198042 5.04943 25.496742\n3 9 0.229013 4.366559 19.066835\n3 13 0.353972 2.825081 7.98108\n3 10 0.257996 3.876026 15.023576\n3 7 0.270829 3.692362 13.633541\n3 0 0.146393 6.830927 46.661562\n3 5 0.127213 7.86081 61.792338\n3 6 0.113836 8.78457 77.168665\n3 1 0.047907 20.87392 435.720551\n3 2 0.028934 34.561005 1194.463096\n3 4 0.05795 17.256279 297.779174\n3 11 0.287624 3.476764 12.087889\n3 12 0.356916 2.801783 7.84999\n3 14 0.34128 2.930144 8.585745\n3 19 0.862919 1.158857 1.34295\n3 18 0.707025 1.414377 2.000461\n3 17 0.612687 1.632154 2.663925\n3 15 0.469555 2.129676 4.535518\n3 16 0.625637 1.59837 2.554787\n0 8 0.335123 2.983978 8.904126\n0 9 0.366094 2.731537 7.461295\n0 13 0.491053 2.036439 4.147086\n0 10 0.395077 2.531151 6.406726\n0 7 0.40778 2.452302 6.013784\n0 3 0.146393 6.830927 46.661562\n0 5 0.264294 3.783661 14.316091\n0 6 0.250787 3.98745 15.89976\n0 1 0.099803 10.019776 100.395905\n0 2 0.118023 8.472919 71.790359\n0 4 0.195031 5.127394 26.290168\n0 11 0.424705 2.354577 5.544034\n0 12 0.493996 2.024306 4.097815\n0 14 0.478361 2.090471 4.370069\n0 19 1 1 1\n0 18 0.844106 1.184685 1.403479\n0 17 0.749768 1.333745 1.778876\n0 15 0.606636 1.648435 2.717338\n0 16 0.762718 1.3111 1.718983\n5 8 0.167174 5.98179 35.78181\n5 9 0.119525 8.366462 69.997687\n5 13 0.246502 4.056755 16.457265\n5 10 0.164657 6.073228 36.884097\n5 7 0.182956 5.465785 29.874806\n5 3 0.127213 7.86081 61.792338\n5 0 0.264294 3.783661 14.316091\n5 6 0.045598 21.931 480.968778\n5 1 0.165808 6.031074 36.373859\n5 2 0.146836 6.810337 46.380692\n5 4 0.069967 14.292411 204.273008\n5 11 0.180154 5.550807 30.811463\n5 12 0.249446 4.008888 16.071183\n5 14 0.23381 4.276971 18.292478\n5 19 0.755449 1.323716 1.752223\n5 18 0.599555 1.667902 2.781898\n5 17 0.505218 1.979345 3.917806\n5 15 0.362085 2.76178 7.62743\n5 16 0.518168 1.929878 3.724428\n6 8 0.199933 5.001679 25.016795\n6 9 0.164249 6.088325 37.067702\n6 13 0.291226 3.433755 11.790674\n6 10 0.209381 4.775982 22.810001\n6 7 0.157261 6.358875 40.435291\n6 3 0.113836 8.78457 77.168665\n6 0 0.250787 3.98745 15.89976\n6 5 0.045598 21.931 480.968778\n6 1 0.1523 6.565968 43.111932\n6 2 0.133328 7.500292 56.254375\n6 4 0.068691 14.558031 211.936253\n6 11 0.224878 4.446857 19.774539\n6 12 0.29417 3.399399 11.555911\n6 14 0.278534 3.590222 12.889695\n6 19 0.800173 1.249729 1.561824\n6 18 0.644279 1.552122 2.409081\n6 17 0.549942 1.818375 3.306487\n6 15 0.406809 2.458155 6.042524\n6 16 0.562891 1.776541 3.156099\n1 8 0.236637 4.225886 17.858116\n1 9 0.267608 3.736811 13.963754\n1 13 0.392567 2.547337 6.488927\n1 10 0.296591 3.371649 11.368015\n1 7 0.309294 3.233172 10.453399\n1 3 0.047907 20.87392 435.720551\n1 0 0.099803 10.019776 100.395905\n1 5 0.165808 6.031074 36.373859\n1 6 0.1523 6.565968 43.111932\n1 2 0.019537 51.185638 2619.96956\n1 4 0.096545 10.357917 107.28645\n1 11 0.326218 3.065432 9.39687\n1 12 0.39551 2.52838 6.392707\n1 14 0.379875 2.632447 6.929776\n1 19 0.901514 1.109246 1.230426\n1 18 0.74562 1.341166 1.798726\n1 17 0.651282 1.535433 2.357555\n1 15 0.50815 1.967924 3.872725\n1 16 0.664232 1.505498 2.266525\n2 8 0.217664 4.594228 21.106932\n2 9 0.248636 4.021951 16.176089\n2 13 0.373594 2.676699 7.164719\n2 10 0.277618 3.602066 12.974877\n2 7 0.290321 3.444457 11.864286\n2 3 0.028934 34.561005 1194.463096\n2 0 0.118023 8.472919 71.790359\n2 5 0.146836 6.810337 46.380692\n2 6 0.133328 7.500292 56.254375\n2 1 0.019537 51.185638 2619.96956\n2 4 0.077572 12.89122 166.183542\n2 11 0.307246 3.254721 10.593207\n2 12 0.376538 2.655776 7.053146\n2 14 0.360902 2.770832 7.677511\n2 19 0.882541 1.133091 1.283896\n2 18 0.726648 1.376183 1.89388\n2 17 0.63231 1.581503 2.501153\n2 15 0.489177 2.044248 4.178951\n2 16 0.64526 1.549764 2.401768\n4 8 0.147367 6.785774 46.046732\n4 9 0.171767 5.821835 33.893763\n4 13 0.297856 3.357331 11.271675\n4 10 0.207321 4.823434 23.265513\n4 7 0.225684 4.430975 19.633541\n4 3 0.05795 17.256279 297.779174\n4 0 0.195031 5.127394 26.290168\n4 5 0.069967 14.292411 204.273008\n4 6 0.068691 14.558031 211.936253\n4 1 0.096545 10.357917 107.28645\n4 2 0.077572 12.89122 166.183542\n4 11 0.231507 4.319521 18.658259\n4 12 0.300799 3.32448 11.052167\n4 14 0.285164 3.50676 12.297363\n4 19 0.806802 1.239461 1.536263\n4 18 0.650909 1.536314 2.36026\n4 17 0.556571 1.796716 3.22819\n4 15 0.413438 2.418739 5.850301\n4 16 0.569521 1.755862 3.083052\n11 8 0.131898 7.581596 57.480603\n11 9 0.103927 9.622145 92.585678\n11 13 0.07454 13.415666 179.980082\n11 10 0.047953 20.853585 434.871996\n11 7 0.255568 3.912846 15.310363\n11 3 0.287624 3.476764 12.087889\n11 0 0.424705 2.354577 5.544034\n11 5 0.180154 5.550807 30.811463\n11 6 0.224878 4.446857 19.774539\n11 1 0.326218 3.065432 9.39687\n11 2 0.307246 3.254721 10.593207\n11 4 0.231507 4.319521 18.658259\n11 12 0.077483 12.906049 166.566106\n11 14 0.06579 15.199795 231.03376\n11 19 0.583487 1.713835 2.937232\n11 18 0.427593 2.338674 5.469394\n11 17 0.333255 3.000705 9.004233\n11 15 0.190123 5.259764 27.665114\n11 16 0.346205 2.888463 8.34322\n12 8 0.196719 5.083381 25.840766\n12 9 0.17682 5.655482 31.984481\n12 13 0.012081 82.771665 6851.148445\n12 10 0.112774 8.867253 78.628178\n12 7 0.328461 3.0445 9.268983\n12 3 0.356916 2.801783 7.84999\n12 0 0.493996 2.024306 4.097815\n12 5 0.249446 4.008888 16.071183\n12 6 0.29417 3.399399 11.555911\n12 1 0.39551 2.52838 6.392707\n12 2 0.376538 2.655776 7.053146\n12 4 0.300799 3.32448 11.052167\n12 11 0.077483 12.906049 166.566106\n12 14 0.04305 23.228547 539.565384\n12 19 0.51449 1.943674 3.777867\n12 18 0.358596 2.788654 7.776593\n12 17 0.264258 3.78418 14.320016\n12 15 0.121126 8.255886 68.159652\n12 16 0.277208 3.6074 13.013335\n14 8 0.1913 5.227401 27.325723\n14 9 0.134911 7.412309 54.94233\n14 13 0.053941 18.538685 343.68286\n14 10 0.107355 9.314917 86.767676\n14 7 0.286552 3.489765 12.178459\n14 3 0.34128 2.930144 8.585745\n14 0 0.478361 2.090471 4.370069\n14 5 0.23381 4.276971 18.292478\n14 6 0.278534 3.590222 12.889695\n14 1 0.379875 2.632447 6.929776\n14 2 0.360902 2.770832 7.677511\n14 4 0.285164 3.50676 12.297363\n14 11 0.06579 15.199795 231.03376\n14 12 0.04305 23.228547 539.565384\n14 19 0.55635 1.797431 3.230759\n14 18 0.400456 2.497155 6.235783\n14 17 0.306118 3.266715 10.671429\n14 15 0.162986 6.135514 37.64453\n14 16 0.319068 3.13413 9.822773\n19 8 0.702723 1.423036 2.025031\n19 9 0.683935 1.462128 2.137818\n19 13 0.509581 1.962397 3.851002\n19 10 0.618778 1.616088 2.611742\n19 7 0.835576 1.196779 1.43228\n19 3 0.862919 1.158857 1.34295\n19 0 1 1 1\n19 5 0.755449 1.323716 1.752223\n19 6 0.800173 1.249729 1.561824\n19 1 0.901514 1.109246 1.230426\n19 2 0.882541 1.133091 1.283896\n19 4 0.806802 1.239461 1.536263\n19 11 0.583487 1.713835 2.937232\n19 12 0.51449 1.943674 3.777867\n19 14 0.55635 1.797431 3.230759\n19 18 0.160045 6.248259 39.040736\n19 17 0.25124 3.980256 15.842441\n19 15 0.478771 2.088682 4.36259\n19 16 0.353693 2.827311 7.993686\n18 8 0.546829 1.828725 3.344233\n18 9 0.528041 1.893793 3.586451\n18 13 0.353687 2.827358 7.993953\n18 10 0.462884 2.160367 4.667187\n18 7 0.679682 1.471275 2.164651\n18 3 0.707025 1.414377 2.000461\n18 0 0.844106 1.184685 1.403479\n18 5 0.599555 1.667902 2.781898\n18 6 0.644279 1.552122 2.409081\n18 1 0.74562 1.341166 1.798726\n18 2 0.726648 1.376183 1.89388\n18 4 0.650909 1.536314 2.36026\n18 11 0.427593 2.338674 5.469394\n18 12 0.358596 2.788654 7.776593\n18 14 0.400456 2.497155 6.235783\n18 19 0.160045 6.248259 39.040736\n18 17 0.095346 10.488083 109.999876\n18 15 0.322877 3.097153 9.592358\n18 16 0.197799 5.055632 25.559411\n17 8 0.452491 2.209987 4.884042\n17 9 0.433703 2.305725 5.316368\n17 13 0.259349 3.855804 14.867228\n17 10 0.368546 2.713362 7.362335\n17 7 0.585345 1.708395 2.918615\n17 3 0.612687 1.632154 2.663925\n17 0 0.749768 1.333745 1.778876\n17 5 0.505218 1.979345 3.917806\n17 6 0.549942 1.818375 3.306487\n17 1 0.651282 1.535433 2.357555\n17 2 0.63231 1.581503 2.501153\n17 4 0.556571 1.796716 3.22819\n17 11 0.333255 3.000705 9.004233\n17 12 0.264258 3.78418 14.320016\n17 14 0.306118 3.266715 10.671429\n17 19 0.25124 3.980256 15.842441\n17 18 0.095346 10.488083 109.999876\n17 15 0.228539 4.375615 19.146003\n17 16 0.149045 6.709402 45.016071\n15 8 0.309359 3.23249 10.448993\n15 9 0.290571 3.441503 11.843945\n15 13 0.116217 8.604601 74.039161\n15 10 0.225414 4.43628 19.680584\n15 7 0.442212 2.261358 5.113739\n15 3 0.469555 2.129676 4.535518\n15 0 0.606636 1.648435 2.717338\n15 5 0.362085 2.76178 7.62743\n15 6 0.406809 2.458155 6.042524\n15 1 0.50815 1.967924 3.872725\n15 2 0.489177 2.044248 4.178951\n15 4 0.413438 2.418739 5.850301\n15 11 0.190123 5.259764 27.665114\n15 12 0.121126 8.255886 68.159652\n15 14 0.162986 6.135514 37.64453\n15 19 0.478771 2.088682 4.36259\n15 18 0.322877 3.097153 9.592358\n15 17 0.228539 4.375615 19.146003\n15 16 0.224747 4.44945 19.797603\n16 8 0.465441 2.148499 4.616047\n16 9 0.446653 2.238875 5.01256\n16 13 0.272299 3.672432 13.486755\n16 10 0.381496 2.621257 6.87099\n16 7 0.598294 1.671418 2.793637\n16 3 0.625637 1.59837 2.554787\n16 0 0.762718 1.3111 1.718983\n16 5 0.518168 1.929878 3.724428\n16 6 0.562891 1.776541 3.156099\n16 1 0.664232 1.505498 2.266525\n16 2 0.64526 1.549764 2.401768\n16 4 0.569521 1.755862 3.083052\n16 11 0.346205 2.888463 8.34322\n16 12 0.277208 3.6074 13.013335\n16 14 0.319068 3.13413 9.822773\n16 19 0.353693 2.827311 7.993686\n16 18 0.197799 5.055632 25.559411\n16 17 0.149045 6.709402 45.016071\n16 15 0.224747 4.44945 19.797603\n
1
0
String

MONITOR
415
750
509
795
alert-activity
alert-level-changes / count locales / ticks
4
1
11

BUTTON
729
750
885
785
toggle-connections
ask connections [set hidden? not hidden?]
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
890
750
1016
784
NIL
toggle-labels
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
457
100
530
134
go-100
repeat 100 [go]
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
1045
55
1218
88
R-clin
R-clin
0
11
3.0
0.01
1
NIL
HORIZONTAL

MONITOR
1125
99
1227
144
clinical-cases
count cases with [clinical?]
0
1
11

MONITOR
1233
100
1315
145
in-hospital
count cases with [in-hospital?]
0
1
11

INPUTBOX
9
472
283
580
alert-levels-control
pessimistic [1 0.8 0.6 0.36]\nrealistic [1 0.72 0.52 0.32]\noptimistic [1 0.64 0.44 0.28]\nother [1 0.8 0.55 0.35]
1
1
String

PLOT
1043
447
1405
741
Daily counts
Days
Count
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"cases" 1.0 0 -2674135 true "" "plot sum [new-cases] of locales"
"infected" 1.0 0 -13345367 true "" "plot sum [new-infected] of locales"
"recovered" 1.0 0 -13840069 true "" "plot sum [new-recovered] of locales"

SWITCH
10
640
152
673
fast-isolation?
fast-isolation?
1
1
-1000

BUTTON
458
139
532
173
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

SLIDER
1229
55
1393
88
p-clinical
p-clinical
0.5
1
0.667
0.001
1
NIL
HORIZONTAL

MONITOR
1328
6
1392
51
R-mean
p-clinical * R-clin + (1 - p-clinical) * 0.5 * R-clin
3
1
11

SLIDER
1124
13
1296
46
new-exposures-arriving
new-exposures-arriving
0
50
0.0
0.01
1
NIL
HORIZONTAL

CHOOSER
10
589
153
634
control-scenario
control-scenario
"optimistic" "realistic" "pessimistic" "other"
1

MONITOR
179
761
261
806
daily tests
sum [item 0 recent-tests] of locales
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
    <enumeratedValueSet variable="log-folder">
      <value value="&quot;staging-area/retest-num-locales&quot;"/>
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
      <value value="false"/>
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

myshape
0.6
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 135 165 150 150
Line -7500403 true 165 165 150 150
@#$#@#$#@
0
@#$#@#$#@
