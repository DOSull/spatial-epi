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


breed [markers marker]

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

  new-exposed
  new-presymptomatic
  new-infected
  new-recovered
  new-dead

  recent-positive-tests   ;; list of recent new positive tests
  recent-tests            ;; list of recent numbers of tests in this locale
  recent-control-levels

  alert-level        ;; local alert level which controls...
  my-alert-indicator ;; local level turtle indicator
  my-control         ;; controls local R0 and
  my-trans-coeff     ;; local transmission rate

  name

  geo-x              ;; geographical coordinates provided
  geo-y              ;; pref. in projected coordinate system
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
  w                  ;; inverse distance weighted _outward_ from each
  gw                 ;; gravity weighted (pop1 x pop2 / dist^2)
  flow-rate-w
  flow-rate-gw
  flow
]


globals [
  n-icu
  cfr-tot

  control-levels         ;; R0 multipliers associated with the alert levels
  flow-levels
  trigger-levels
  scripted-events   ;; list of timestamp alert-level changes

  mean-R0           ;; pop weighted mean R0
  mean-trans-coeff  ;; pop weighted mean trans coeff

  total-exposed            ;; exposed
  total-presymptomatic     ;; infected but not symptomatic
  total-infected           ;; infected and symptomatic
  total-recovered          ;; recovered
  total-dead               ;; dead

  total-new-exposed
  total-new-presymptomatic
  total-new-infected
  total-new-recovered
  total-new-dead

  log-header-file-name
  log-file-name
  full-file-name
  date-time
  model-name
  labels-on?

  alert-level-changes

  start-time ;; ticks at which setup is complete
]

to setup
  clear-all
  reset-ticks

  setup-visualization-defaults

  ;; random seed for PRNG for replicability
  if use-seed? [random-seed seed]

  set control-levels initialise-control-levels
  set flow-levels read-from-string alert-levels-flow
  set trigger-levels read-from-string alert-level-triggers
  set scripted-events map [x -> map [y -> read-from-string y] split-string x " "] split-string script "\n"

  setup-locales
  setup-levels

  ;; some rescaling of sizes for visibility
  let size-adjust 0.25 / sqrt (count locales / count patches)
  ask (turtle-set locales levels) [ set size size * size-adjust ]
  ask connections [ set thickness thickness * size-adjust ]

  ;; initial exposures
  setup-cases
  enact-alert-levels

  update-statistics

  if not netlogo-web? and log-all-locales? [
    initialise-logging
  ]

  set labels-on? true
  redraw
  paint-land green 4

  set start-time ticks
end

;; main initialisation of locales
to setup-locales
  ;; using data in various strings at end of the model code
  if setup-method = "Random landscape" [
    setup-locales-parametrically
    setup-connections-parametrically
    stop
  ]
  if position "NZ DHBs" setup-method = 0 [
    setup-locales-from-string get-locales-data "DHBs"
    setup-connections-from-string get-connections-data true
    stop
  ]
  if position "NZ TAs" setup-method = 0 [
    setup-locales-from-string get-locales-data "TAs"
    setup-connections-from-string get-connections-data false
    stop
  ]
  setup-locales-from-string get-locales-data "costa-rica"
  setup-connections-parametrically
end

;; initialisation of cases
to setup-cases
  ;; from NZ MoH data at DHB level
  ifelse member? "MoH data" setup-method [
    initialise-cases-from-string get-cases-data
  ]
  [ ;; otherwise random initialisation
    ;; this will be uniform by population
    ifelse initialise-by-burn-in? [
      burn-in-model
    ]
    [
      add-arrivals initial-infected
    ]
  ]
end


to burn-in-model
  let save-initial-alert-level initial-alert-level
  let save-alert-policy alert-policy
  let save-new-exposures-arriving new-exposures-arriving
  let save-pop-test-rate pop-test-rate
  let save-inf-test-rate inf-test-rate
  set-to-unprepared

  carefully [
    while [total-burden-remaining < initial-infected] [
      run-one-day true
      tick
    ]
    set-parameters save-initial-alert-level save-alert-policy save-new-exposures-arriving save-pop-test-rate save-inf-test-rate
  ]
  [
    set-parameters save-initial-alert-level save-alert-policy save-new-exposures-arriving save-pop-test-rate save-inf-test-rate
  ]
end


to set-to-unprepared
  set-parameters 1 "static" 5 0 0
end


to set-parameters [a-level policy arrivals pop-test-r inf-test-r]
  ask locales [
    set alert-level a-level
  ]
  enact-alert-levels
  set alert-policy policy
  set new-exposures-arriving arrivals
  set pop-test-rate pop-test-r
  set inf-test-rate inf-test-r
end




;; add new infections in relative proportions
;; dictated by the relative rate of movement between states
to initial-infect
  let rel-frequencies (list (1 / exp-to-presymp) (1 / presymp-to-inf) (1 / inf-to-rec))
  let cum-rel-frequencies cumulative-sum rel-frequencies
  let total-rel-frequencies last cum-rel-frequencies
  let i random total-rel-frequencies
  let choice length filter [x -> x < i] cum-rel-frequencies
  set susceptible susceptible - 1
  if choice = 0 [
    set exposed exposed + 1
    stop
  ]
  if choice = 1 [
    set presymptomatic presymptomatic + 1
    stop
  ]
  set infected infected + 1
end




to-report cumulative-sum [lst]
  let starter sublist lst 0 1
  report reduce [ [a b] -> lput (last a + b) a] (fput starter but-first lst)
end

to update-statistics
  set mean-R0 sum [my-control * R0 * susceptible] of locales / sum [susceptible] of locales
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
  if (total-burden-remaining = 0 and ticks > 0) [
    update-plots
    tick
    stop
  ]
  reset-timer
  ;; run one day of the model
  run-one-day false
  if timer > 5 [
    output-print "Epidemic appears out of control"
    output-print "or model has slowed for some other"
    output-print "reason, stopping execution"
    stop
  ]

  if not netlogo-web? and log-all-locales? [
    file-open full-file-name
    ask locales [
      file-print output-locale
    ]
    file-close
  ]

  tick
end


to run-one-day [burn-in?]
  ;; add some new cases if appropriate
  add-arrivals random-poisson new-exposures-arriving

  ask locales [
    set-control-level
    calculate-flows
    spread
  ]

  update-statistics
  update-testing

  if not burn-in? [
    ifelse alert-policy = "scripted" [
      let elapsed-time ticks - start-time
      let events filter [x -> item 0 x = elapsed-time] scripted-events
      if debug? [ show events ]
      if length events > 0 [
        set initial-alert-level last first events
        change-alert-levels
      ]
    ]
    [
      ;; change alert levels
      if ticks >= (start-lifting-quarantine + start-time) and
      ticks >= time-horizon and
      (ticks - start-lifting-quarantine - start-time) mod time-horizon = 0 [
        change-alert-levels
      ]
    ]
  ]
  redraw
end

to set-control-level
  ifelse response-time > 0 [
    set recent-control-levels fput (item (alert-level - 1) control-levels) recent-control-levels
    set my-control mean sublist recent-control-levels 0 min (list response-time length recent-control-levels)
  ]
  [
    set my-control item (alert-level - 1) control-levels
  ]
  set my-trans-coeff get-transmission-coeff (my-control * R0)
end

;; susceptible weighted infection so
;; that higher population locales are more
;; likely to receive exposures
to add-arrivals [n]
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


;; if this reporter falls to 0 then the model may as well stop
;; and we've beaten this bastarding virus!
to-report total-burden-remaining
  report total-infected + total-presymptomatic + total-exposed
end



;; ------------------------------
;; ALERT LEVELS
;; ------------------------------
;; In all cases alert levels are 'sticky' downwards, i.e. they won't jump immediately
;; to the lowest possible level consistent positive tesing rates
to change-alert-levels
  if alert-policy = "static" [ ;; there will be no changes unless user has changed the initial-alert-level slider
    ask locales [
      let new-level initial-alert-level
      if new-level != alert-level [
        set alert-level new-level
        set alert-level-changes alert-level-changes + 1
      ]
    ]
    enact-alert-levels
    stop
  ]

  if alert-policy = "local-random" [ ;; move alert levels up or down a level at random
    ask locales [
      let alert-level-change one-of (list -1 0 1)
      let new-level clamp (alert-level + alert-level-change) 1 4
      if new-level != alert-level [
        set alert-level new-level
        set alert-level-changes alert-level-changes + 1
      ]
    ]
    enact-alert-levels
    stop
  ]


  if alert-policy = "local" [
    ask locales [ ;; change alert levels locally based on local positive test results
      let positive-test-rate get-positive-test-rate
      let new-level get-new-alert-level positive-test-rate
      if new-level != alert-level [
        ifelse new-level < alert-level [ ;; only go down one level
          set alert-level clamp (alert-level - 1) 1 4
        ]
        [
          set alert-level new-level
        ]
        set alert-level-changes alert-level-changes + 1
      ]
    ]
    enact-alert-levels
    stop
  ]

  if alert-policy = "global-max" [ ;; base global changes on the maximum local positive test rate
    let positive-test-rate max [get-positive-test-rate] of locales
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
    enact-alert-levels
    stop
  ]

  if alert-policy = "global-mean" [ ;; base global changes on population weighted (ie total population) based rates
    let positive-test-rate sum [get-positive-test-rate * pop-0] of locales / sum [pop-0] of locales
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
    enact-alert-levels
    stop
  ]
end

to update-testing
  ask locales [
    if debug? [
      show (list ticks (length recent-tests) (length recent-positive-tests))
    ]
    let new-infected-tests random-poisson (infected * inf-test-rate)
    let new-tests new-infected-tests + random-poisson (susceptible * pop-test-rate)
    set recent-tests fput new-tests recent-tests

    let positive-tests random-binomial new-infected-tests (1 - false-negative-rate)
    set recent-positive-tests fput positive-tests recent-positive-tests
  ]
end


to-report get-positive-test-rate
  let current-tests recent-total recent-tests
  let current-positive-tests recent-total recent-positive-tests
  ifelse current-tests > 0 [
    report current-positive-tests / current-tests
  ]
  [ ;; if there have been no tests, then caution dictates a high rate
    report 1
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

to enact-alert-levels
  ask locales [
;    set mandated-control item (alert-level - 1) control-levels
;    ifelse response-time > 0 and length my-recent-control-levels > 0 [
;      let available-data min (list (response-time + 1) length my-recent-control-levels)
;      let recent-mean mean sublist my-recent-control-levels 0 available-data
;      set my-trans-coeff get-transmission-coeff (recent-mean * R0)
;    ]
;    [
;      set my-trans-coeff get-transmission-coeff (my-control * R0)
;    ]
  ]
  ask levels [
    set alert-level [alert-level] of my-locale
    draw-level
  ]
  ask connections [
    set flow-rate-w w * min [item (alert-level - 1) flow-levels] of both-ends
    set flow-rate-gw gw * min [item (alert-level - 1) flow-levels] of both-ends
  ]
end


;to set-transmission-parameters [a]
;  set alert-level a
;  set my-R0 item alert-level R0-levels
;  set my-trans-coeff get-transmission-coeff my-R0
;end

to-report get-transmission-coeff [R]
  report R / (rel-inf-presymp / presymp-to-inf + 1 / inf-to-rec)
end



to spread
  let immunity-correction 1 / (pop-0 - dead)
  set new-exposed random-binomial (susceptible) (get-effective-infected * my-trans-coeff * immunity-correction) +
                  random-binomial (susceptible) (get-effective-presymptomatic * my-trans-coeff * rel-inf-presymp * immunity-correction)
  set new-presymptomatic random-binomial exposed exp-to-presymp
  set new-infected random-binomial presymptomatic presymp-to-inf

  let no-longer-infected random-binomial infected inf-to-rec
  set new-recovered random-binomial no-longer-infected (1 - cfr-tot)
  set new-dead no-longer-infected - new-recovered

;  show (list new-exposed new-presymptomatic new-infected no-longer-infected new-recovered new-dead)

  ;; update all the stocks
  expose new-exposed
  presym new-presymptomatic
  infect new-infected
  recover new-recovered
  kill new-dead
end

to expose [n]
  ;show word "expose " n
  set susceptible susceptible - n
  set exposed exposed + n
end

to presym [n]
  ;show word "presym " n
  set exposed exposed - n
  set presymptomatic presymptomatic + n
end

to infect [n]
  ;show word "infect " n
  set presymptomatic presymptomatic - n
  set infected infected + n
end

to recover [n]
  ;show word "recover " n
  set infected infected - n
  set recovered recovered + n
end

to kill [n]
  ;show word "kill " n
  set infected infected - n
  set pop-0 pop-0 - n
  set dead dead + n
end

to calculate-flows
  ask my-out-connections [
    ifelse gravity-weight? [
      set flow random-poisson ([presymptomatic] of myself * flow-rate-gw)
    ]
    [
      set flow random-poisson ([presymptomatic] of myself * flow-rate-w)
    ]
  ]
end

to-report get-effective-infected
  report infected
;  report infected + flow-rate * sum [my-flow-rate * w * [infected] of other-end] of my-in-connections
end

to-report get-effective-presymptomatic
  report presymptomatic - sum [flow] of my-out-connections + sum [flow] of my-in-connections
  ; report presymptomatic + flow-rate * sum [my-flow-rate * w * [presymptomatic] of other-end] of my-in-connections
end

;; --------------------------------------------------------------
;; NOTE
;; using random-poisson approximation for efficiency when n large
;; --------------------------------------------------------------
;;
;to-report random-binomial [n p]
;  if p = 0 [ report 0 ]
;  if p = 1 [ report n ]
;  if n * p >= 10 and p < 0.25 [
;    report random-poisson (n * p)
;  ]
;  report length filter [x -> x < p] n-values n [random-float 1]
;end


;; Based on code from https://stackoverflow.com/questions/23561551/a-efficient-binomial-random-number-generator-code-in-java#23574723
;; which implements one of Devroye's algorithms - and not the super complicated BTPE one of Kach... something
to-report random-binomial [n p]
; the Java code
;public static int getBinomial(int n, double p) {
;   double log_q = Math.log(1.0 - p);
;   int x = 0;
;   double sum = 0;
;   for(;;) {
;      sum += Math.log(Math.random()) / (n - x);
;      if(sum < log_q) {
;         return x;
;      }
;      x++;
;   }
  ; need to trap p = 0 and p = 1
  if p = 1 [ report n ]
  if p = 0 [ report 0 ]
  let ln-q ln (1 - p)
  let x 0
  let s 0
  ; also need to avoid x = n
  while [x < n] [
    set s s + ln (random-float 1) / (n - x)
    if s < ln-q [
      report x
    ]
    set x x + 1
  ]
  report x
end


;; --------------------------------
;; initialisation stuff
;; --------------------------------
;; locales
to setup-locales-parametrically
  let pop-mean population / num-locales
  let pop-var pop-sd-multiplier * pop-sd-multiplier * pop-mean * pop-mean

  let alpha (pop-mean * pop-mean / pop-var)
  let lambda (pop-mean / pop-var)

  create-locales num-locales [
    let xy random-xy-with-buffer 0.1
    setxy item 0 xy item 1 xy
    set geo-x xcor
    set geo-y ycor
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
    set exposed 0
    set presymptomatic 0
    set infected 0
    set recovered 0
    set dead 0

    set new-exposed 0
    set new-presymptomatic 0
    set new-infected 0
    set new-recovered 0
    set new-dead 0

    set recent-tests []
    set recent-positive-tests []
    set recent-control-levels n-values response-time [x -> 1]

    set alert-level initial-alert-level
    set my-control item (alert-level - 1) control-levels
    set my-trans-coeff get-transmission-coeff (my-control * R0)
  ]
end


;; reads from a multiline string where each line is "id name x y pop\n"
;; lines are assumed to be sorted from id 0 by
to setup-locales-from-string [s]
  ;; make it into a list of strings one per locale
  ;; first line is a header for humans to read, but we don't need it
  let locales-data but-first split-string s "\n"

  ;; we need to do some preprocessing of the x y coordinates for rescaling purposes
  let xs map [x -> read-from-string item 2 split-string x " "] locales-data
  let ys map [y -> read-from-string item 3 split-string y " "] locales-data
  let bbox get-bbox xs ys
  create-locales length locales-data

  (foreach locales-data sort locales [ [line loc] ->
    let parameters split-string line " "
    ask loc [
      set name item 1 parameters
      set label name
      set label-color black
      set pop-0 read-from-string item 4 parameters
      set geo-x read-from-string item 2 parameters
      set geo-y read-from-string item 3 parameters
      let xy rescale-xy (list geo-x geo-y) bbox 1
      setxy item 0 xy item 1 xy
;      setxy (rescale geo-x min-x max-x 0.1 true)
;            (rescale geo-y min-y max-y 0.1 false)
    ]
  ])
  initialise-locales
end


to-report get-bbox [x-coords y-coords]
  let min-x min x-coords
  let max-x max x-coords
  let min-y min y-coords
  let max-y max y-coords

  report (list min-x min-y max-x max-y)
end


;; rescale supplied xy coordinate pair (a list [x y]) to specified bbox
;; with a buffer to the Netlogo world boundary
to-report rescale-xy [xy bbox buffer]
  let old-x-range abs item 2 bbox - item 0 bbox
  let old-y-range abs item 3 bbox - item 1 bbox
  let old-x-min min (list item 0 bbox item 2 bbox)
  let old-y-min min (list item 1 bbox item 3 bbox)
  let new-x-range max-pxcor - min-pxcor - (2 * buffer)
  let new-y-range max-pycor - min-pycor - (2 * buffer)
  let scale min (list (new-x-range / old-x-range) (new-y-range / old-y-range))
  set new-x-range old-x-range * scale
  set new-y-range old-y-range * scale
  let new-x-min 0 - new-x-range / 2
  let new-y-min 0 - new-y-range / 2
  report (list (new-x-min + (item 0 xy - old-x-min) * scale)
               (new-y-min + (item 1 xy - old-y-min) * scale))
end


;; rescale input x or y coordinate to the world coordinate system
;; z is the value to be rescaled, z-min and z-max the range of the z values
;; buffer is a proportion of the range inset from the world edge to keep the
;; rescaled coordinats
;; x? is true if we are dealing with x coordinates, false if y coordinates
to-report rescale [z min-z max-z buffer x?]
  let new-range ifelse-value x? [(max-pxcor - min-pxcor) * (1 - buffer)] [(max-pycor - min-pycor) * (1 - buffer)]
  let new-min ifelse-value x? [min-pxcor + buffer / 2 * new-range] [min-pycor + buffer / 2 * new-range]
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

;; initialise from a multiline string "timestamp DHB\n"
to initialise-cases-from-string [s]
  let cases-data but-first split-string s "\n"

  foreach cases-data [ c ->
    let parameters split-string c " "
    let dhb-name item 1 parameters
    ask one-of locales with [name = dhb-name] [
      initial-infect
    ]
  ]
end


;; their connections
to setup-connections-parametrically
  ask locales [
    create-pair-of-connections self nearest-non-neighbour
  ]
  let num-links (8 * count locales)
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



to setup-connections-from-string [s]
  let edges but-first split-string s "\n"
  foreach edges [ edge ->
    let parameters split-string edge " "
    let v1 locale (read-from-string item 0 parameters)
    let v2 locale (read-from-string item 1 parameters)
    let weight read-from-string item 2 parameters
    ask v1 [
      create-connection-to v2 [
        initialise-connection 1 / weight
      ]
    ]
  ]
  let max-d max list (max-connection-distance * 1000) (max [min [geo-distance] of my-out-connections] of locales)
  ask connections with [geo-distance > max-d] [ die ]
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
  set color [0 115 230 96]
  set w 1 / d
  let gd geo-distance
  set gw 1 / gd / gd
  set thickness w
end

to reweight-connections
  ask locales [
    let p1 pop-0
    ask my-out-connections [
      set gw gw * p1 * [pop-0] of other-end
    ]
    let total-w sum [w] of my-out-connections
    let total-gw sum [gw] of my-out-connections
    let w-correction 1 / total-w
    let gw-correction 1 / total-gw
    let thickness-correction 1 / total-w
    ask my-out-connections [
      set w w * w-correction
      set gw gw * gw-correction
      set thickness thickness * thickness-correction
    ]
  ]
end


;; connection length
to-report geo-distance
  let ends sort both-ends
  let xs map [v -> [geo-x] of v] ends
  let ys map [v -> [geo-y] of v] ends
  report sqrt ((item 0 xs - item 1 xs) ^ 2 + (item 0 ys - item 1 ys) ^ 2)
end



;; ----------------------------
;; Drawing stuff
;; ----------------------------
to setup-visualization-defaults
  ;; background colour i.e. sea
  ask patches [set pcolor cyan + 1]

  ;; shapes of the different turtle types
  set-default-shape locales "circle"
  set-default-shape levels "square 3"
  set-default-shape connections "myshape"
  set-default-shape markers "circle 3"
end

to redraw
  clear-drawing
  ask locales [
    draw-locale
  ]
  ask levels [
    draw-level
  ]
end

to draw-locale
  set color scale-color red recovered pop-0 0
  if (exposed + presymptomatic + infected) > 0 [
    hatch 1 [
      set breed markers
      set size size * 1.1
      set color black
      stamp
      die
    ]
  ]
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
  ;;
  let lake patches with [shade-of? pcolor cyan and count neighbors4 with [pcolor = c] >= 3]
  ask lake [
    set pcolor c
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
  set model-name "distributed-seir-RC-logging"
  set date-time replace date-and-time "." ":"
  let base-file-name (word model-name "-" date-time)
  set base-file-name join-list split-string base-file-name " " "-"
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
  report join-list (list ticks who name pop-0 susceptible exposed presymptomatic infected recovered dead new-exposed new-presymptomatic new-infected new-recovered new-dead (first recent-tests) alert-level) ","
end

to-report output-locale-header
  report "ticks,who,name,pop.0,susceptible,exposed,presymptomatic,infected,recovered,dead,new.exposed,new.presymptomatic,new.infected,new.recovered,new.dead,new.tests,alert.level"
end

to-report log-file-header
  let parameters (list "name,value")
  set parameters lput join-list (list "model.name" model-name) "," parameters
  set parameters lput join-list (list "log.file.name" log-file-name) "," parameters
  set parameters lput join-list (list "log.all.locales?" log-all-locales?) "," parameters
  set parameters lput join-list (list "log.folder" log-folder) "," parameters

  set parameters lput join-list (list "use.seed?" use-seed?) "," parameters
  set parameters lput join-list (list "seed" seed) "," parameters
  set parameters lput join-list (list "setup.method" setup-method) "," parameters
  set parameters lput join-list (list "gravity.weight?" gravity-weight?) "," parameters
  set parameters lput join-list (list "max-connection.distance" max-connection-distance) "," parameters

  set parameters lput join-list (list "initial.infected" initial-infected) "," parameters
  set parameters lput join-list (list "initial.alert.level" initial-alert-level) "," parameters
  set parameters lput join-list (list "alert.policy" alert-policy) "," parameters
  set parameters lput join-list (list "start.lifting.quarantine" start-lifting-quarantine) "," parameters
  set parameters lput join-list (list "time.horizon" time-horizon) "," parameters

  set parameters lput join-list (list "alert.levels.control" alert-levels-control) "," parameters
  set parameters lput join-list (list "alert.levels.flow" alert-levels-flow) "," parameters

  set parameters lput join-list (list "pop.test.rate" pop-test-rate) "," parameters
  set parameters lput join-list (list "inf.test.rate" inf-test-rate) "," parameters
  set parameters lput join-list (list "false.negative.rate" false-negative-rate) "," parameters

  set parameters lput join-list (list "population" population) "," parameters
  set parameters lput join-list (list "num.locales" num-locales) "," parameters
  set parameters lput join-list (list "pop.sd.multiplier" pop-sd-multiplier) "," parameters

  set parameters lput join-list (list "exp.to.presymp" exp-to-presymp) "," parameters
  set parameters lput join-list (list "presymp.to.inf" presymp-to-inf) "," parameters
  set parameters lput join-list (list "rel.inf.presymp" rel-inf-presymp) "," parameters
  set parameters lput join-list (list "inf.to.rec" inf-to-rec) "," parameters
  set parameters lput join-list (list "p.hosp" p-hosp) "," parameters
  set parameters lput join-list (list "p.icu" p-icu) "," parameters
  set parameters lput join-list (list "icu.cap" icu-cap) "," parameters
  set parameters lput join-list (list "cfr.0" cfr-0) "," parameters
  set parameters lput join-list (list "cfr.1" cfr-1) "," parameters

  set parameters lput join-list (list "min.pxcor" min-pxcor) "," parameters
  set parameters lput join-list (list "max.pxcor" max-pxcor) "," parameters
  set parameters lput join-list (list "min.pycor" min-pycor) "," parameters
  set parameters lput join-list (list "max.pycor" max-pycor) "," parameters

  report join-list parameters "\n"
end




to-report get-locales-data [dataset]
  if dataset = "DHBs" [
    report join-list [ "ID name x y pop"
      "0 Northland 1700348 6061149 178500"
      "1 Waitemata 1749476 5947203 622350"
      "2 Auckland 1757594 5920685 538840"
      "3 Counties.Manukau 1766973 5885442 557790"
      "4 Waikato 1810656 5816377 417130"
      "5 Lakes 1878016 5754014 110040"
      "6 Bay.of.Plenty 1896923 5815743 236820"
      "7 Tairawhiti 2038766 5713797 48965"
      "8 Taranaki 1699580 5662304 119640"
      "9 Hawkes.Bay 1934136 5612568 165360"
      "10 Whanganui 1784458 5579703 64565"
      "11 MidCentral 1816860 5524911 178240"
      "12 Hutt.Valley 1763986 5437753 149270"
      "13 Capital.and.Coast 1753070 5437211 316620"
      "14 Wairarapa 1817576 5457737 44845"
      "15 Nelson.Marlborough 1636138 5423667 150280"
      "16 West.Coast 1458495 5311071 32445"
      "17 Canterbury 1566860 5181157 562870"
      "18 South.Canterbury 1455463 5085777 60025"
      "19 Southern 1343425 4917902 328230"
    ] "\n"
  ]
  if dataset = "TAs" [
    report join-list [ "ID name x y pop"
      "35 Porirua.City 1754803 5444679 59100"
      "36 Upper.Hutt.City 1773987 5445388 46000"
      "37 Lower.Hutt.City 1760219 5434877 108700"
      "38 Wellington.City 1749100 5428463 210400"
      "39 Masterton.District 1824544 5463335 26800"
      "40 Carterton.District 1812507 5455416 9690"
      "41 South.Wairarapa.District 1800827 5443083 11100"
      "42 Tasman.District 1605286 5434837 54800"
      "0 Far.North.District 1661974 6104072 68500"
      "1 Whangarei.District 1721923 6044508 96000"
      "2 Kaipara.District 1696226 6013981 24100"
      "3 Thames-Coromandel.District 1837677 5900170 31500"
      "4 Hauraki.District 1839725 5862162 21000"
      "5 Waikato.District 1785683 5838225 79900"
      "6 Matamata-Piako.District 1835194 5825452 36000"
      "7 Hamilton.City 1800871 5815452 169500"
      "8 Waipa.District 1810124 5797080 56200"
      "9 Ōtorohanga.District 1790241 5774620 10500"
      "10 South.Waikato.District 1849099 5771898 25100"
      "23 Central.Hawke's.Bay.District 1904903 5568810 14850"
      "24 New.Plymouth.District 1695862 5676313 84400"
      "25 Stratford.District 1711134 5645373 9860"
      "26 South.Taranaki.District 1708671 5619311 28600"
      "27 Ruapehu.District 1802853 5669639 12750"
      "28 Whanganui.District 1775205 5577837 47300"
      "29 Rangitikei.District 1809324 5568325 15750"
      "30 Manawatu.District 1815966 5543552 31700"
      "31 Palmerston.North.City 1822162 5529590 88300"
      "32 Tararua.District 1854975 5533899 18650"
      "33 Horowhenua.District 1793597 5504781 35000"
      "34 Kapiti.Coast.District 1771493 5472047 56000"
      "11 Waitomo.District 1785591 5753009 9490"
      "12 Taupo.District 1861966 5710310 39300"
      "13 Western.Bay.of.Plenty.District 1879557 5826857 53900"
      "14 Tauranga.City 1881487 5825032 144700"
      "15 Rotorua.District 1885518 5774443 75100"
      "16 Whakatane.District 1946637 5786008 37100"
      "17 Kawerau.District 1924725 5778172 7490"
      "18 Ōpōtiki.District 1979259 5787601 9720"
      "19 Gisborne.District 2038766 5713797 49300"
      "20 Wairoa.District 1986371 5670970 8680"
      "21 Hastings.District 1930906 5605005 85000"
      "22 Napier.City 1936851 5621688 65000"
      "43 Nelson.City 1621899 5428781 52900"
      "44 Marlborough.District 1679129 5408151 49200"
      "45 Kaikoura.District 1655830 5305721 4110"
      "46 Buller.District 1490283 5372488 9840"
      "47 Grey.District 1456087 5299058 13750"
      "48 Westland.District 1419121 5250185 8960"
      "49 Hurunui.District 1587134 5244241 12950"
      "50 Waimakariri.District 1568013 5202021 62800"
      "51 Christchurch.City 1570759 5179833 385500"
      "52 Selwyn.District 1544814 5171648 65600"
      "53 Ashburton.District 1499844 5141411 34800"
      "54 Timaru.District 1460651 5088516 47900"
      "55 Mackenzie.District 1392391 5111453 5140"
      "56 Waimate.District 1446643 5045173 8080"
      "57 Waitaki.District 1435380 5003721 23200"
      "58 Central.Otago.District 1315931 4989515 23100"
      "59 Queenstown-Lakes.District 1271815 5018598 41700"
      "60 Dunedin.City 1405074 4917548 131700"
      "61 Clutha.District 1350001 4881249 18350"
      "62 Southland.District 1227217 4891543 32100"
      "63 Gore.District 1285924 4885546 12800"
      "64 Invercargill.City 1242280 4848618 56200"
      "65 Auckland 1757446 5921056 1642800" ] "\n"
  ]
  if dataset = "costa-rica" [
    report join-list [ "ID canton x y pop2020"
      "0 Abangares 394965 1136868 20016"
      "1 Acosta 482021 1083629 21976"
      "2 Quepos 482269 1042921 33069"
      "3 Alajuela 476248 1107595 314209"
      "4 Alajuelita 488825 1095147 94548"
      "5 Alvarado 520955 1096595 15433"
      "6 Aserrí 489646 1091029 63529"
      "7 Atenas 457976 1103921 29340"
      "8 Bagaces 362664 1164390 24130"
      "9 Barva 486510 1108050 47002"
      "10 Belén 479626 1103384 26459"
      "11 Buenos.Aires 573256 1011776 53436"
      "12 Cañas 380566 1153107 32685"
      "13 Carrillo 330292 1155315 45939"
      "14 Cartago 509139 1090999 164121"
      "15 Corredores 616183 956498 52419"
      "16 Coto.Brus 613669 975401 44308"
      "17 Curridabat 495735 1096405 79577"
      "18 Desamparados 493103 1095059 245208"
      "19 Dota 503444 1067312 7948"
      "20 El.Guarco 507311 1089155 46304"
      "21 Escazú 484407 1097022 70054"
      "22 Esparza 426754 1105080 38183"
      "23 Flores 483403 1106116 24886"
      "24 Garabito 431097 1063412 26028"
      "25 Goicoechea 493999 1099969 138525"
      "26 Golfito 591347 955437 45573"
      "27 Grecia 465900 1113747 93845"
      "28 Guácimo 534693 1128802 55128"
      "29 Guatuso 410655 1179594 19236"
      "30 Heredia 487209 1105347 143208"
      "31 Hojancha 344675 1112469 7998"
      "32 Jiménez 527843 1094633 16321"
      "33 La.Cruz 321536 1224210 27090"
      "34 La.Unión 501644 1095514 112508"
      "35 León.Cortés.Castro 495610 1070814 13769"
      "36 Liberia 343164 1176152 76969"
      "37 Limón 605984 1104055 99836"
      "38 Los.Chiles 421691 1220124 33689"
      "39 Matina 577698 1114245 46379"
      "40 Montes.de.Oca 494155 1098310 62533"
      "41 Montes.de.Oro 419904 1116252 14323"
      "42 Mora 473560 1096261 30318"
      "43 Moravia 494518 1102058 62669"
      "44 Nandayure 362769 1105389 11787"
      "45 Naranjo 457988 1116828 48803"
      "46 Nicoya 341096 1122688 56591"
      "47 Oreamuno 510605 1091276 49972"
      "48 Orotina 442581 1095989 23786"
      "49 Osa 552414 990515 31139"
      "50 Palmares 452504 1111305 143117"
      "51 Paraíso 514807 1087898 40928"
      "52 Parrita 463899 1052832 62941"
      "53 Pérez.Zeledón 532369 1036697 20199"
      "54 Poás 473238 1114202 34006"
      "55 Pococí 523740 1129718 150664"
      "56 Puntarenas 408572 1103248 140102"
      "57 Puriscal 465265 1089171 37983"
      "58 San.Carlos 452819 1141834 200151"
      "59 San.Isidro 493970 1107649 23230"
      "60 San.José 491354 1098295 347398"
      "61 San.Mateo 442859 1098723 7141"
      "62 San.Pablo 489493 1105346 31200"
      "63 San.Rafael 489189 1107374 55269"
      "64 San.Ramón 448486 1115401 93872"
      "65 Santa.Ana 479897 1098376 60453"
      "66 Santa.Barbara 482644 1109865 42778"
      "67 Santa.Cruz 326601 1134967 68939"
      "68 Santo.Domingo 490132 1103502 49045"
      "69 Sarapiquí 499391 1156127 83015"
      "70 Siquirres 554037 1116476 64923"
      "71 Talamanca 626069 1064605 43153"
      "72 Tarrazú 497682 1068121 18535"
      "73 Tibás 491162 1101659 84873"
      "74 Tilarán 394175 1157518 21749"
      "75 Turrialba 534717 1095382 73659"
      "76 Turrubares 451749 1095515 6871"
      "77 Upala 388919 1205285 54055"
      "78 Sarchí 461679 1115547 71663"
      "79 Vázquez.de.Coronado 499185 1103057 22166"
      "80 Zarcero 457151 1126234 14341" ] "\n"
  ]
end

to-report get-connections-data [dhbs-or-tas]
  ifelse dhbs-or-tas [
    report join-list [ "ID1 ID2 cost weight weight2"
      "8 9 0.227069 4.40394 19.394684"
      "8 13 0.193776 5.160594 26.631733"
      "8 10 0.093423 10.704056 114.576815"
      "8 7 0.345958 2.890528 8.355151"
      "8 3 0.198042 5.04943 25.496742"
      "8 0 0.335123 2.983978 8.904126"
      "8 5 0.167174 5.98179 35.78181"
      "8 6 0.199933 5.001679 25.016795"
      "8 1 0.236637 4.225886 17.858116"
      "8 2 0.217664 4.594228 21.106932"
      "8 4 0.147367 6.785774 46.046732"
      "8 11 0.131898 7.581596 57.480603"
      "8 12 0.196719 5.083381 25.840766"
      "8 14 0.1913 5.227401 27.325723"
      "8 19 0.702723 1.423036 2.025031"
      "8 18 0.546829 1.828725 3.344233"
      "8 17 0.452491 2.209987 4.884042"
      "8 15 0.309359 3.23249 10.448993"
      "8 16 0.465441 2.148499 4.616047"
      "9 8 0.227069 4.40394 19.394684"
      "9 13 0.174988 5.714684 32.657614"
      "9 10 0.143124 6.986926 48.817133"
      "9 7 0.154221 6.484195 42.044779"
      "9 3 0.229013 4.366559 19.066835"
      "9 0 0.366094 2.731537 7.461295"
      "9 5 0.119525 8.366462 69.997687"
      "9 6 0.164249 6.088325 37.067702"
      "9 1 0.267608 3.736811 13.963754"
      "9 2 0.248636 4.021951 16.176089"
      "9 4 0.171767 5.821835 33.893763"
      "9 11 0.103927 9.622145 92.585678"
      "9 12 0.17682 5.655482 31.984481"
      "9 14 0.134911 7.412309 54.94233"
      "9 19 0.683935 1.462128 2.137818"
      "9 18 0.528041 1.893793 3.586451"
      "9 17 0.433703 2.305725 5.316368"
      "9 15 0.290571 3.441503 11.843945"
      "9 16 0.446653 2.238875 5.01256"
      "13 8 0.193776 5.160594 26.631733"
      "13 9 0.174988 5.714684 32.657614"
      "13 10 0.109831 9.104883 82.898897"
      "13 7 0.326629 3.061574 9.373238"
      "13 3 0.353972 2.825081 7.98108"
      "13 0 0.491053 2.036439 4.147086"
      "13 5 0.246502 4.056755 16.457265"
      "13 6 0.291226 3.433755 11.790674"
      "13 1 0.392567 2.547337 6.488927"
      "13 2 0.373594 2.676699 7.164719"
      "13 4 0.297856 3.357331 11.271675"
      "13 11 0.07454 13.415666 179.980082"
      "13 12 0.012081 82.771665 6851.148445"
      "13 14 0.053941 18.538685 343.68286"
      "13 19 0.509581 1.962397 3.851002"
      "13 18 0.353687 2.827358 7.993953"
      "13 17 0.259349 3.855804 14.867228"
      "13 15 0.116217 8.604601 74.039161"
      "13 16 0.272299 3.672432 13.486755"
      "10 8 0.093423 10.704056 114.576815"
      "10 9 0.143124 6.986926 48.817133"
      "10 13 0.109831 9.104883 82.898897"
      "10 7 0.294766 3.392522 11.509203"
      "10 3 0.257996 3.876026 15.023576"
      "10 0 0.395077 2.531151 6.406726"
      "10 5 0.164657 6.073228 36.884097"
      "10 6 0.209381 4.775982 22.810001"
      "10 1 0.296591 3.371649 11.368015"
      "10 2 0.277618 3.602066 12.974877"
      "10 4 0.207321 4.823434 23.265513"
      "10 11 0.047953 20.853585 434.871996"
      "10 12 0.112774 8.867253 78.628178"
      "10 14 0.107355 9.314917 86.767676"
      "10 19 0.618778 1.616088 2.611742"
      "10 18 0.462884 2.160367 4.667187"
      "10 17 0.368546 2.713362 7.362335"
      "10 15 0.225414 4.43628 19.680584"
      "10 16 0.381496 2.621257 6.87099"
      "7 8 0.345958 2.890528 8.355151"
      "7 9 0.154221 6.484195 42.044779"
      "7 13 0.326629 3.061574 9.373238"
      "7 10 0.294766 3.392522 11.509203"
      "7 3 0.270829 3.692362 13.633541"
      "7 0 0.40778 2.452302 6.013784"
      "7 5 0.182956 5.465785 29.874806"
      "7 6 0.157261 6.358875 40.435291"
      "7 1 0.309294 3.233172 10.453399"
      "7 2 0.290321 3.444457 11.864286"
      "7 4 0.225684 4.430975 19.633541"
      "7 11 0.255568 3.912846 15.310363"
      "7 12 0.328461 3.0445 9.268983"
      "7 14 0.286552 3.489765 12.178459"
      "7 19 0.835576 1.196779 1.43228"
      "7 18 0.679682 1.471275 2.164651"
      "7 17 0.585345 1.708395 2.918615"
      "7 15 0.442212 2.261358 5.113739"
      "7 16 0.598294 1.671418 2.793637"
      "3 8 0.198042 5.04943 25.496742"
      "3 9 0.229013 4.366559 19.066835"
      "3 13 0.353972 2.825081 7.98108"
      "3 10 0.257996 3.876026 15.023576"
      "3 7 0.270829 3.692362 13.633541"
      "3 0 0.146393 6.830927 46.661562"
      "3 5 0.127213 7.86081 61.792338"
      "3 6 0.113836 8.78457 77.168665"
      "3 1 0.047907 20.87392 435.720551"
      "3 2 0.028934 34.561005 1194.463096"
      "3 4 0.05795 17.256279 297.779174"
      "3 11 0.287624 3.476764 12.087889"
      "3 12 0.356916 2.801783 7.84999"
      "3 14 0.34128 2.930144 8.585745"
      "3 19 0.862919 1.158857 1.34295"
      "3 18 0.707025 1.414377 2.000461"
      "3 17 0.612687 1.632154 2.663925"
      "3 15 0.469555 2.129676 4.535518"
      "3 16 0.625637 1.59837 2.554787"
      "0 8 0.335123 2.983978 8.904126"
      "0 9 0.366094 2.731537 7.461295"
      "0 13 0.491053 2.036439 4.147086"
      "0 10 0.395077 2.531151 6.406726"
      "0 7 0.40778 2.452302 6.013784"
      "0 3 0.146393 6.830927 46.661562"
      "0 5 0.264294 3.783661 14.316091"
      "0 6 0.250787 3.98745 15.89976"
      "0 1 0.099803 10.019776 100.395905"
      "0 2 0.118023 8.472919 71.790359"
      "0 4 0.195031 5.127394 26.290168"
      "0 11 0.424705 2.354577 5.544034"
      "0 12 0.493996 2.024306 4.097815"
      "0 14 0.478361 2.090471 4.370069"
      "0 19 1 1 1"
      "0 18 0.844106 1.184685 1.403479"
      "0 17 0.749768 1.333745 1.778876"
      "0 15 0.606636 1.648435 2.717338"
      "0 16 0.762718 1.3111 1.718983"
      "5 8 0.167174 5.98179 35.78181"
      "5 9 0.119525 8.366462 69.997687"
      "5 13 0.246502 4.056755 16.457265"
      "5 10 0.164657 6.073228 36.884097"
      "5 7 0.182956 5.465785 29.874806"
      "5 3 0.127213 7.86081 61.792338"
      "5 0 0.264294 3.783661 14.316091"
      "5 6 0.045598 21.931 480.968778"
      "5 1 0.165808 6.031074 36.373859"
      "5 2 0.146836 6.810337 46.380692"
      "5 4 0.069967 14.292411 204.273008"
      "5 11 0.180154 5.550807 30.811463"
      "5 12 0.249446 4.008888 16.071183"
      "5 14 0.23381 4.276971 18.292478"
      "5 19 0.755449 1.323716 1.752223"
      "5 18 0.599555 1.667902 2.781898"
      "5 17 0.505218 1.979345 3.917806"
      "5 15 0.362085 2.76178 7.62743"
      "5 16 0.518168 1.929878 3.724428"
      "6 8 0.199933 5.001679 25.016795"
      "6 9 0.164249 6.088325 37.067702"
      "6 13 0.291226 3.433755 11.790674"
      "6 10 0.209381 4.775982 22.810001"
      "6 7 0.157261 6.358875 40.435291"
      "6 3 0.113836 8.78457 77.168665"
      "6 0 0.250787 3.98745 15.89976"
      "6 5 0.045598 21.931 480.968778"
      "6 1 0.1523 6.565968 43.111932"
      "6 2 0.133328 7.500292 56.254375"
      "6 4 0.068691 14.558031 211.936253"
      "6 11 0.224878 4.446857 19.774539"
      "6 12 0.29417 3.399399 11.555911"
      "6 14 0.278534 3.590222 12.889695"
      "6 19 0.800173 1.249729 1.561824"
      "6 18 0.644279 1.552122 2.409081"
      "6 17 0.549942 1.818375 3.306487"
      "6 15 0.406809 2.458155 6.042524"
      "6 16 0.562891 1.776541 3.156099"
      "1 8 0.236637 4.225886 17.858116"
      "1 9 0.267608 3.736811 13.963754"
      "1 13 0.392567 2.547337 6.488927"
      "1 10 0.296591 3.371649 11.368015"
      "1 7 0.309294 3.233172 10.453399"
      "1 3 0.047907 20.87392 435.720551"
      "1 0 0.099803 10.019776 100.395905"
      "1 5 0.165808 6.031074 36.373859"
      "1 6 0.1523 6.565968 43.111932"
      "1 2 0.019537 51.185638 2619.96956"
      "1 4 0.096545 10.357917 107.28645"
      "1 11 0.326218 3.065432 9.39687"
      "1 12 0.39551 2.52838 6.392707"
      "1 14 0.379875 2.632447 6.929776"
      "1 19 0.901514 1.109246 1.230426"
      "1 18 0.74562 1.341166 1.798726"
      "1 17 0.651282 1.535433 2.357555"
      "1 15 0.50815 1.967924 3.872725"
      "1 16 0.664232 1.505498 2.266525"
      "2 8 0.217664 4.594228 21.106932"
      "2 9 0.248636 4.021951 16.176089"
      "2 13 0.373594 2.676699 7.164719"
      "2 10 0.277618 3.602066 12.974877"
      "2 7 0.290321 3.444457 11.864286"
      "2 3 0.028934 34.561005 1194.463096"
      "2 0 0.118023 8.472919 71.790359"
      "2 5 0.146836 6.810337 46.380692"
      "2 6 0.133328 7.500292 56.254375"
      "2 1 0.019537 51.185638 2619.96956"
      "2 4 0.077572 12.89122 166.183542"
      "2 11 0.307246 3.254721 10.593207"
      "2 12 0.376538 2.655776 7.053146"
      "2 14 0.360902 2.770832 7.677511"
      "2 19 0.882541 1.133091 1.283896"
      "2 18 0.726648 1.376183 1.89388"
      "2 17 0.63231 1.581503 2.501153"
      "2 15 0.489177 2.044248 4.178951"
      "2 16 0.64526 1.549764 2.401768"
      "4 8 0.147367 6.785774 46.046732"
      "4 9 0.171767 5.821835 33.893763"
      "4 13 0.297856 3.357331 11.271675"
      "4 10 0.207321 4.823434 23.265513"
      "4 7 0.225684 4.430975 19.633541"
      "4 3 0.05795 17.256279 297.779174"
      "4 0 0.195031 5.127394 26.290168"
      "4 5 0.069967 14.292411 204.273008"
      "4 6 0.068691 14.558031 211.936253"
      "4 1 0.096545 10.357917 107.28645"
      "4 2 0.077572 12.89122 166.183542"
      "4 11 0.231507 4.319521 18.658259"
      "4 12 0.300799 3.32448 11.052167"
      "4 14 0.285164 3.50676 12.297363"
      "4 19 0.806802 1.239461 1.536263"
      "4 18 0.650909 1.536314 2.36026"
      "4 17 0.556571 1.796716 3.22819"
      "4 15 0.413438 2.418739 5.850301"
      "4 16 0.569521 1.755862 3.083052"
      "11 8 0.131898 7.581596 57.480603"
      "11 9 0.103927 9.622145 92.585678"
      "11 13 0.07454 13.415666 179.980082"
      "11 10 0.047953 20.853585 434.871996"
      "11 7 0.255568 3.912846 15.310363"
      "11 3 0.287624 3.476764 12.087889"
      "11 0 0.424705 2.354577 5.544034"
      "11 5 0.180154 5.550807 30.811463"
      "11 6 0.224878 4.446857 19.774539"
      "11 1 0.326218 3.065432 9.39687"
      "11 2 0.307246 3.254721 10.593207"
      "11 4 0.231507 4.319521 18.658259"
      "11 12 0.077483 12.906049 166.566106"
      "11 14 0.06579 15.199795 231.03376"
      "11 19 0.583487 1.713835 2.937232"
      "11 18 0.427593 2.338674 5.469394"
      "11 17 0.333255 3.000705 9.004233"
      "11 15 0.190123 5.259764 27.665114"
      "11 16 0.346205 2.888463 8.34322"
      "12 8 0.196719 5.083381 25.840766"
      "12 9 0.17682 5.655482 31.984481"
      "12 13 0.012081 82.771665 6851.148445"
      "12 10 0.112774 8.867253 78.628178"
      "12 7 0.328461 3.0445 9.268983"
      "12 3 0.356916 2.801783 7.84999"
      "12 0 0.493996 2.024306 4.097815"
      "12 5 0.249446 4.008888 16.071183"
      "12 6 0.29417 3.399399 11.555911"
      "12 1 0.39551 2.52838 6.392707"
      "12 2 0.376538 2.655776 7.053146"
      "12 4 0.300799 3.32448 11.052167"
      "12 11 0.077483 12.906049 166.566106"
      "12 14 0.04305 23.228547 539.565384"
      "12 19 0.51449 1.943674 3.777867"
      "12 18 0.358596 2.788654 7.776593"
      "12 17 0.264258 3.78418 14.320016"
      "12 15 0.121126 8.255886 68.159652"
      "12 16 0.277208 3.6074 13.013335"
      "14 8 0.1913 5.227401 27.325723"
      "14 9 0.134911 7.412309 54.94233"
      "14 13 0.053941 18.538685 343.68286"
      "14 10 0.107355 9.314917 86.767676"
      "14 7 0.286552 3.489765 12.178459"
      "14 3 0.34128 2.930144 8.585745"
      "14 0 0.478361 2.090471 4.370069"
      "14 5 0.23381 4.276971 18.292478"
      "14 6 0.278534 3.590222 12.889695"
      "14 1 0.379875 2.632447 6.929776"
      "14 2 0.360902 2.770832 7.677511"
      "14 4 0.285164 3.50676 12.297363"
      "14 11 0.06579 15.199795 231.03376"
      "14 12 0.04305 23.228547 539.565384"
      "14 19 0.55635 1.797431 3.230759"
      "14 18 0.400456 2.497155 6.235783"
      "14 17 0.306118 3.266715 10.671429"
      "14 15 0.162986 6.135514 37.64453"
      "14 16 0.319068 3.13413 9.822773"
      "19 8 0.702723 1.423036 2.025031"
      "19 9 0.683935 1.462128 2.137818"
      "19 13 0.509581 1.962397 3.851002"
      "19 10 0.618778 1.616088 2.611742"
      "19 7 0.835576 1.196779 1.43228"
      "19 3 0.862919 1.158857 1.34295"
      "19 0 1 1 1"
      "19 5 0.755449 1.323716 1.752223"
      "19 6 0.800173 1.249729 1.561824"
      "19 1 0.901514 1.109246 1.230426"
      "19 2 0.882541 1.133091 1.283896"
      "19 4 0.806802 1.239461 1.536263"
      "19 11 0.583487 1.713835 2.937232"
      "19 12 0.51449 1.943674 3.777867"
      "19 14 0.55635 1.797431 3.230759"
      "19 18 0.160045 6.248259 39.040736"
      "19 17 0.25124 3.980256 15.842441"
      "19 15 0.478771 2.088682 4.36259"
      "19 16 0.353693 2.827311 7.993686"
      "18 8 0.546829 1.828725 3.344233"
      "18 9 0.528041 1.893793 3.586451"
      "18 13 0.353687 2.827358 7.993953"
      "18 10 0.462884 2.160367 4.667187"
      "18 7 0.679682 1.471275 2.164651"
      "18 3 0.707025 1.414377 2.000461"
      "18 0 0.844106 1.184685 1.403479"
      "18 5 0.599555 1.667902 2.781898"
      "18 6 0.644279 1.552122 2.409081"
      "18 1 0.74562 1.341166 1.798726"
      "18 2 0.726648 1.376183 1.89388"
      "18 4 0.650909 1.536314 2.36026"
      "18 11 0.427593 2.338674 5.469394"
      "18 12 0.358596 2.788654 7.776593"
      "18 14 0.400456 2.497155 6.235783"
      "18 19 0.160045 6.248259 39.040736"
      "18 17 0.095346 10.488083 109.999876"
      "18 15 0.322877 3.097153 9.592358"
      "18 16 0.197799 5.055632 25.559411"
      "17 8 0.452491 2.209987 4.884042"
      "17 9 0.433703 2.305725 5.316368"
      "17 13 0.259349 3.855804 14.867228"
      "17 10 0.368546 2.713362 7.362335"
      "17 7 0.585345 1.708395 2.918615"
      "17 3 0.612687 1.632154 2.663925"
      "17 0 0.749768 1.333745 1.778876"
      "17 5 0.505218 1.979345 3.917806"
      "17 6 0.549942 1.818375 3.306487"
      "17 1 0.651282 1.535433 2.357555"
      "17 2 0.63231 1.581503 2.501153"
      "17 4 0.556571 1.796716 3.22819"
      "17 11 0.333255 3.000705 9.004233"
      "17 12 0.264258 3.78418 14.320016"
      "17 14 0.306118 3.266715 10.671429"
      "17 19 0.25124 3.980256 15.842441"
      "17 18 0.095346 10.488083 109.999876"
      "17 15 0.228539 4.375615 19.146003"
      "17 16 0.149045 6.709402 45.016071"
      "15 8 0.309359 3.23249 10.448993"
      "15 9 0.290571 3.441503 11.843945"
      "15 13 0.116217 8.604601 74.039161"
      "15 10 0.225414 4.43628 19.680584"
      "15 7 0.442212 2.261358 5.113739"
      "15 3 0.469555 2.129676 4.535518"
      "15 0 0.606636 1.648435 2.717338"
      "15 5 0.362085 2.76178 7.62743"
      "15 6 0.406809 2.458155 6.042524"
      "15 1 0.50815 1.967924 3.872725"
      "15 2 0.489177 2.044248 4.178951"
      "15 4 0.413438 2.418739 5.850301"
      "15 11 0.190123 5.259764 27.665114"
      "15 12 0.121126 8.255886 68.159652"
      "15 14 0.162986 6.135514 37.64453"
      "15 19 0.478771 2.088682 4.36259"
      "15 18 0.322877 3.097153 9.592358"
      "15 17 0.228539 4.375615 19.146003"
      "15 16 0.224747 4.44945 19.797603"
      "16 8 0.465441 2.148499 4.616047"
      "16 9 0.446653 2.238875 5.01256"
      "16 13 0.272299 3.672432 13.486755"
      "16 10 0.381496 2.621257 6.87099"
      "16 7 0.598294 1.671418 2.793637"
      "16 3 0.625637 1.59837 2.554787"
      "16 0 0.762718 1.3111 1.718983"
      "16 5 0.518168 1.929878 3.724428"
      "16 6 0.562891 1.776541 3.156099"
      "16 1 0.664232 1.505498 2.266525"
      "16 2 0.64526 1.549764 2.401768"
      "16 4 0.569521 1.755862 3.083052"
      "16 11 0.346205 2.888463 8.34322"
      "16 12 0.277208 3.6074 13.013335"
      "16 14 0.319068 3.13413 9.822773"
      "16 19 0.353693 2.827311 7.993686"
      "16 18 0.197799 5.055632 25.559411"
      "16 17 0.149045 6.709402 45.016071"
      "16 15 0.224747 4.44945 19.797603" ] "\n"
  ]
  [
    report join-list [ "ID1 ID2 cost"
      "0 1 0.5"
      "0 2 0.5"
      "0 65 1"
      "1 0 0.5"
      "1 2 0.5"
      "1 65 1"
      "2 0 0.5"
      "2 1 0.5"
      "2 65 0.5"
      "2 4 1"
      "2 5 1"
      "3 4 0.5"
      "3 5 1"
      "3 6 1"
      "3 13 1"
      "3 65 1"
      "4 3 0.5"
      "4 5 0.5"
      "4 6 0.5"
      "4 13 0.5"
      "4 65 0.5"
      "4 7 1"
      "4 8 1"
      "4 9 1"
      "4 10 1"
      "4 14 1"
      "4 15 1"
      "4 16 1"
      "4 2 1"
      "5 4 0.5"
      "5 6 0.5"
      "5 7 0.5"
      "5 8 0.5"
      "5 9 0.5"
      "5 65 0.5"
      "5 3 1"
      "5 13 1"
      "5 10 1"
      "5 11 1"
      "5 12 1"
      "5 2 1"
      "6 4 0.5"
      "6 5 0.5"
      "6 8 0.5"
      "6 10 0.5"
      "6 13 0.5"
      "6 3 1"
      "6 65 1"
      "6 7 1"
      "6 9 1"
      "6 12 1"
      "6 15 1"
      "6 14 1"
      "6 16 1"
      "7 5 0.5"
      "7 8 0.5"
      "7 4 1"
      "7 6 1"
      "7 9 1"
      "7 65 1"
      "7 10 1"
      "8 5 0.5"
      "8 6 0.5"
      "8 7 0.5"
      "8 9 0.5"
      "8 10 0.5"
      "8 4 1"
      "8 65 1"
      "8 13 1"
      "8 11 1"
      "8 12 1"
      "8 15 1"
      "9 5 0.5"
      "9 8 0.5"
      "9 10 0.5"
      "9 11 0.5"
      "9 12 0.5"
      "9 4 1"
      "9 6 1"
      "9 7 1"
      "9 65 1"
      "9 13 1"
      "9 15 1"
      "9 24 1"
      "9 27 1"
      "9 29 1"
      "9 16 1"
      "9 20 1"
      "9 21 1"
      "10 6 0.5"
      "10 8 0.5"
      "10 9 0.5"
      "10 12 0.5"
      "10 13 0.5"
      "10 15 0.5"
      "10 4 1"
      "10 5 1"
      "10 7 1"
      "10 11 1"
      "10 27 1"
      "10 29 1"
      "10 16 1"
      "10 20 1"
      "10 21 1"
      "10 14 1"
      "11 9 0.5"
      "11 24 0.5"
      "11 27 0.5"
      "11 12 0.5"
      "11 5 1"
      "11 8 1"
      "11 10 1"
      "11 25 1"
      "11 26 1"
      "11 28 1"
      "11 29 1"
      "11 15 1"
      "11 16 1"
      "11 20 1"
      "11 21 1"
      "12 9 0.5"
      "12 10 0.5"
      "12 27 0.5"
      "12 29 0.5"
      "12 11 0.5"
      "12 15 0.5"
      "12 16 0.5"
      "12 20 0.5"
      "12 21 0.5"
      "12 5 1"
      "12 8 1"
      "12 6 1"
      "12 13 1"
      "12 24 1"
      "12 25 1"
      "12 28 1"
      "12 23 1"
      "12 30 1"
      "12 17 1"
      "12 18 1"
      "12 19 1"
      "12 22 1"
      "13 4 0.5"
      "13 6 0.5"
      "13 10 0.5"
      "13 14 0.5"
      "13 15 0.5"
      "13 16 0.5"
      "13 3 1"
      "13 5 1"
      "13 65 1"
      "13 8 1"
      "13 9 1"
      "13 12 1"
      "13 17 1"
      "13 18 1"
      "13 19 1"
      "13 20 1"
      "13 21 1"
      "14 13 0.5"
      "14 4 1"
      "14 6 1"
      "14 10 1"
      "14 15 1"
      "14 16 1"
      "15 10 0.5"
      "15 12 0.5"
      "15 13 0.5"
      "15 16 0.5"
      "15 6 1"
      "15 8 1"
      "15 9 1"
      "15 27 1"
      "15 29 1"
      "15 11 1"
      "15 20 1"
      "15 21 1"
      "15 4 1"
      "15 14 1"
      "15 17 1"
      "15 18 1"
      "15 19 1"
      "16 12 0.5"
      "16 13 0.5"
      "16 15 0.5"
      "16 17 0.5"
      "16 18 0.5"
      "16 19 0.5"
      "16 20 0.5"
      "16 21 0.5"
      "16 9 1"
      "16 10 1"
      "16 27 1"
      "16 29 1"
      "16 11 1"
      "16 4 1"
      "16 6 1"
      "16 14 1"
      "16 23 1"
      "16 22 1"
      "17 16 0.5"
      "17 12 1"
      "17 13 1"
      "17 15 1"
      "17 18 1"
      "17 19 1"
      "17 20 1"
      "17 21 1"
      "18 16 0.5"
      "18 19 0.5"
      "18 12 1"
      "18 13 1"
      "18 15 1"
      "18 17 1"
      "18 20 1"
      "18 21 1"
      "19 16 0.5"
      "19 18 0.5"
      "19 20 0.5"
      "19 12 1"
      "19 13 1"
      "19 15 1"
      "19 17 1"
      "19 21 1"
      "20 12 0.5"
      "20 16 0.5"
      "20 19 0.5"
      "20 21 0.5"
      "20 9 1"
      "20 10 1"
      "20 27 1"
      "20 29 1"
      "20 11 1"
      "20 15 1"
      "20 13 1"
      "20 17 1"
      "20 18 1"
      "20 23 1"
      "20 22 1"
      "21 23 0.5"
      "21 29 0.5"
      "21 12 0.5"
      "21 16 0.5"
      "21 20 0.5"
      "21 22 0.5"
      "21 30 1"
      "21 32 1"
      "21 27 1"
      "21 28 1"
      "21 9 1"
      "21 10 1"
      "21 11 1"
      "21 15 1"
      "21 13 1"
      "21 17 1"
      "21 18 1"
      "21 19 1"
      "22 21 0.5"
      "22 23 1"
      "22 29 1"
      "22 12 1"
      "22 16 1"
      "22 20 1"
      "23 29 0.5"
      "23 30 0.5"
      "23 32 0.5"
      "23 21 0.5"
      "23 27 1"
      "23 28 1"
      "23 12 1"
      "23 31 1"
      "23 33 1"
      "23 39 1"
      "23 16 1"
      "23 20 1"
      "23 22 1"
      "24 25 0.5"
      "24 26 0.5"
      "24 27 0.5"
      "24 11 0.5"
      "24 28 1"
      "24 29 1"
      "24 12 1"
      "24 9 1"
      "25 24 0.5"
      "25 26 0.5"
      "25 27 0.5"
      "25 28 0.5"
      "25 11 1"
      "25 29 1"
      "25 12 1"
      "26 24 0.5"
      "26 25 0.5"
      "26 28 0.5"
      "26 27 1"
      "26 11 1"
      "26 29 1"
      "27 24 0.5"
      "27 25 0.5"
      "27 28 0.5"
      "27 29 0.5"
      "27 11 0.5"
      "27 12 0.5"
      "27 26 1"
      "27 23 1"
      "27 30 1"
      "27 21 1"
      "27 9 1"
      "27 10 1"
      "27 15 1"
      "27 16 1"
      "27 20 1"
      "28 25 0.5"
      "28 26 0.5"
      "28 27 0.5"
      "28 29 0.5"
      "28 24 1"
      "28 11 1"
      "28 12 1"
      "28 23 1"
      "28 30 1"
      "28 21 1"
      "29 23 0.5"
      "29 27 0.5"
      "29 28 0.5"
      "29 30 0.5"
      "29 12 0.5"
      "29 21 0.5"
      "29 32 1"
      "29 24 1"
      "29 25 1"
      "29 11 1"
      "29 26 1"
      "29 31 1"
      "29 33 1"
      "29 9 1"
      "29 10 1"
      "29 15 1"
      "29 16 1"
      "29 20 1"
      "29 22 1"
      "30 23 0.5"
      "30 29 0.5"
      "30 31 0.5"
      "30 32 0.5"
      "30 33 0.5"
      "30 21 1"
      "30 27 1"
      "30 28 1"
      "30 12 1"
      "30 39 1"
      "30 40 1"
      "30 34 1"
      "31 30 0.5"
      "31 32 0.5"
      "31 33 0.5"
      "31 23 1"
      "31 29 1"
      "31 39 1"
      "31 40 1"
      "31 34 1"
      "32 39 0.5"
      "32 23 0.5"
      "32 30 0.5"
      "32 31 0.5"
      "32 33 0.5"
      "32 40 1"
      "32 29 1"
      "32 21 1"
      "32 34 1"
      "33 39 0.5"
      "33 40 0.5"
      "33 30 0.5"
      "33 31 0.5"
      "33 32 0.5"
      "33 34 0.5"
      "33 41 1"
      "33 23 1"
      "33 29 1"
      "33 35 1"
      "33 36 1"
      "34 35 0.5"
      "34 36 0.5"
      "34 40 0.5"
      "34 41 0.5"
      "34 33 0.5"
      "34 37 1"
      "34 38 1"
      "34 39 1"
      "34 30 1"
      "34 31 1"
      "34 32 1"
      "35 36 0.5"
      "35 37 0.5"
      "35 38 0.5"
      "35 34 0.5"
      "35 41 1"
      "35 44 1"
      "35 40 1"
      "35 33 1"
      "36 35 0.5"
      "36 37 0.5"
      "36 41 0.5"
      "36 34 0.5"
      "36 38 1"
      "36 40 1"
      "36 33 1"
      "37 35 0.5"
      "37 36 0.5"
      "37 38 0.5"
      "37 41 0.5"
      "37 34 1"
      "37 44 1"
      "37 40 1"
      "38 35 0.5"
      "38 37 0.5"
      "38 44 0.5"
      "38 36 1"
      "38 34 1"
      "38 41 1"
      "38 42 1"
      "38 43 1"
      "38 45 1"
      "38 49 1"
      "39 40 0.5"
      "39 32 0.5"
      "39 33 0.5"
      "39 41 1"
      "39 34 1"
      "39 23 1"
      "39 30 1"
      "39 31 1"
      "40 39 0.5"
      "40 41 0.5"
      "40 33 0.5"
      "40 34 0.5"
      "40 32 1"
      "40 36 1"
      "40 37 1"
      "40 30 1"
      "40 31 1"
      "40 35 1"
      "41 36 0.5"
      "41 37 0.5"
      "41 40 0.5"
      "41 34 0.5"
      "41 35 1"
      "41 38 1"
      "41 39 1"
      "41 33 1"
      "42 43 0.5"
      "42 44 0.5"
      "42 46 0.5"
      "42 49 0.5"
      "42 45 1"
      "42 38 1"
      "42 47 1"
      "42 48 1"
      "42 50 1"
      "42 52 1"
      "43 42 0.5"
      "43 44 0.5"
      "43 46 1"
      "43 49 1"
      "43 45 1"
      "43 38 1"
      "44 42 0.5"
      "44 43 0.5"
      "44 45 0.5"
      "44 49 0.5"
      "44 38 0.5"
      "44 46 1"
      "44 47 1"
      "44 48 1"
      "44 50 1"
      "44 52 1"
      "44 35 1"
      "44 37 1"
      "45 44 0.5"
      "45 49 0.5"
      "45 42 1"
      "45 43 1"
      "45 38 1"
      "45 46 1"
      "45 47 1"
      "45 48 1"
      "45 50 1"
      "45 52 1"
      "46 42 0.5"
      "46 47 0.5"
      "46 49 0.5"
      "46 43 1"
      "46 44 1"
      "46 48 1"
      "46 45 1"
      "46 50 1"
      "46 52 1"
      "47 46 0.5"
      "47 48 0.5"
      "47 49 0.5"
      "47 42 1"
      "47 52 1"
      "47 53 1"
      "47 55 1"
      "47 57 1"
      "47 59 1"
      "47 62 1"
      "47 44 1"
      "47 45 1"
      "47 50 1"
      "48 47 0.5"
      "48 49 0.5"
      "48 52 0.5"
      "48 53 0.5"
      "48 55 0.5"
      "48 57 0.5"
      "48 59 0.5"
      "48 62 0.5"
      "48 46 1"
      "48 42 1"
      "48 44 1"
      "48 45 1"
      "48 50 1"
      "48 51 1"
      "48 54 1"
      "48 56 1"
      "48 58 1"
      "48 60 1"
      "48 61 1"
      "48 63 1"
      "48 64 1"
      "49 42 0.5"
      "49 44 0.5"
      "49 45 0.5"
      "49 46 0.5"
      "49 47 0.5"
      "49 48 0.5"
      "49 50 0.5"
      "49 52 0.5"
      "49 43 1"
      "49 38 1"
      "49 53 1"
      "49 55 1"
      "49 57 1"
      "49 59 1"
      "49 62 1"
      "49 51 1"
      "50 49 0.5"
      "50 51 0.5"
      "50 52 0.5"
      "50 42 1"
      "50 44 1"
      "50 45 1"
      "50 46 1"
      "50 47 1"
      "50 48 1"
      "50 53 1"
      "51 50 0.5"
      "51 52 0.5"
      "51 49 1"
      "51 48 1"
      "51 53 1"
      "52 48 0.5"
      "52 49 0.5"
      "52 50 0.5"
      "52 51 0.5"
      "52 53 0.5"
      "52 47 1"
      "52 55 1"
      "52 57 1"
      "52 59 1"
      "52 62 1"
      "52 42 1"
      "52 44 1"
      "52 45 1"
      "52 46 1"
      "52 54 1"
      "53 48 0.5"
      "53 52 0.5"
      "53 54 0.5"
      "53 55 0.5"
      "53 47 1"
      "53 49 1"
      "53 57 1"
      "53 59 1"
      "53 62 1"
      "53 50 1"
      "53 51 1"
      "53 56 1"
      "54 53 0.5"
      "54 55 0.5"
      "54 56 0.5"
      "54 48 1"
      "54 52 1"
      "54 57 1"
      "55 48 0.5"
      "55 53 0.5"
      "55 54 0.5"
      "55 56 0.5"
      "55 57 0.5"
      "55 47 1"
      "55 49 1"
      "55 52 1"
      "55 59 1"
      "55 62 1"
      "55 58 1"
      "55 60 1"
      "56 54 0.5"
      "56 55 0.5"
      "56 57 0.5"
      "56 53 1"
      "56 48 1"
      "56 58 1"
      "56 59 1"
      "56 60 1"
      "57 48 0.5"
      "57 55 0.5"
      "57 56 0.5"
      "57 58 0.5"
      "57 59 0.5"
      "57 60 0.5"
      "57 47 1"
      "57 49 1"
      "57 52 1"
      "57 53 1"
      "57 62 1"
      "57 54 1"
      "57 61 1"
      "58 57 0.5"
      "58 59 0.5"
      "58 60 0.5"
      "58 61 0.5"
      "58 62 0.5"
      "58 48 1"
      "58 55 1"
      "58 56 1"
      "58 63 1"
      "58 64 1"
      "59 48 0.5"
      "59 57 0.5"
      "59 58 0.5"
      "59 62 0.5"
      "59 47 1"
      "59 49 1"
      "59 52 1"
      "59 53 1"
      "59 55 1"
      "59 56 1"
      "59 60 1"
      "59 61 1"
      "59 63 1"
      "59 64 1"
      "60 57 0.5"
      "60 58 0.5"
      "60 61 0.5"
      "60 48 1"
      "60 55 1"
      "60 56 1"
      "60 59 1"
      "60 62 1"
      "60 63 1"
      "61 58 0.5"
      "61 60 0.5"
      "61 62 0.5"
      "61 63 0.5"
      "61 57 1"
      "61 59 1"
      "61 48 1"
      "61 64 1"
      "62 48 0.5"
      "62 58 0.5"
      "62 59 0.5"
      "62 61 0.5"
      "62 63 0.5"
      "62 64 0.5"
      "62 47 1"
      "62 49 1"
      "62 52 1"
      "62 53 1"
      "62 55 1"
      "62 57 1"
      "62 60 1"
      "63 61 0.5"
      "63 62 0.5"
      "63 58 1"
      "63 60 1"
      "63 48 1"
      "63 59 1"
      "63 64 1"
      "64 62 0.5"
      "64 48 1"
      "64 58 1"
      "64 59 1"
      "64 61 1"
      "64 63 1"
      "65 2 0.5"
      "65 4 0.5"
      "65 5 0.5"
      "65 0 1"
      "65 1 1"
      "65 3 1"
      "65 6 1"
      "65 13 1"
      "65 7 1"
      "65 8 1"
      "65 9 1"
    ] "\n"
  ]
end

to-report get-cases-data
  report join-list [ "t DHB"
    "0 Southern"
    "0 Waitemata"
    "-1 Canterbury"
    "-1 Counties.Manukau"
    "-1 Counties.Manukau"
    "-1 Counties.Manukau"
    "-1 Counties.Manukau"
    "-1 South.Canterbury"
    "-2 Canterbury"
    "-2 Northland"
    "-2 Tairawhiti"
    "-2 Waitemata"
    "-3 Canterbury"
    "-3 Canterbury"
    "-3 Capital.and.Coast"
    "-3 Hawkes.Bay"
    "-3 Hawkes.Bay"
    "-3 Lakes"
    "-3 Lakes"
    "-3 Waitemata"
    "-3 Waitemata"
    "-3 Waitemata"
    "-3 Waitemata"
    "-4 Auckland"
    "-4 Counties.Manukau"
    "-4 Northland"
    "-4 South.Canterbury"
    "-4 Waikato"
    "-4 Waikato"
    "-4 Waikato"
    "-4 Waitemata"
    "-4 Waitemata"
    "-4 Waitemata"
    "-4 Waitemata"
    "-4 Waitemata"
    "-5 Auckland"
    "-5 Bay.of.Plenty"
    "-5 Counties.Manukau"
    "-5 Northland"
    "-5 Waitemata"
    "-5 Waitemata"
    "-5 Waitemata"
    "-5 Waitemata"
    "-6 Auckland"
    "-6 Auckland"
    "-6 Auckland"
    "-6 Auckland"
    "-6 Auckland"
    "-6 Auckland"
    "-6 Auckland"
    "-6 Auckland"
    "-6 Canterbury"
    "-6 Canterbury"
    "-6 Canterbury"
    "-6 Canterbury"
    "-6 Canterbury"
    "-6 Capital.and.Coast"
    "-6 Capital.and.Coast"
    "-6 Counties.Manukau"
    "-6 Counties.Manukau"
    "-6 Counties.Manukau"
    "-6 Counties.Manukau"
    "-6 Counties.Manukau"
    "-6 Hawkes.Bay"
    "-6 Southern"
    "-6 Southern"
    "-6 Waikato"
    "-6 Waikato"
    "-6 Waikato"
    "-6 Waitemata"
    "-6 Waitemata"
    "-6 Waitemata"
    "-6 Waitemata"
    "-7 Auckland"
    "-7 Auckland"
    "-7 Auckland"
    "-7 Canterbury"
    "-7 Capital.and.Coast"
    "-7 Capital.and.Coast"
    "-7 Capital.and.Coast"
    "-7 Counties.Manukau"
    "-7 Counties.Manukau"
    "-7 Hawkes.Bay"
    "-7 South.Canterbury"
    "-7 Southern"
    "-7 Southern"
    "-7 Waitemata"
    "-7 Waitemata"
    "-7 Waitemata"
    "-7 Waitemata"
    "-8 Auckland"
    "-8 Bay.of.Plenty"
    "-8 Canterbury"
    "-8 Counties.Manukau"
    "-8 Hawkes.Bay"
    "-8 Northland"
    "-8 Northland"
    "-8 Southern"
    "-8 Southern"
    "-8 Southern"
    "-8 Southern"
    "-8 Southern"
    "-8 Southern"
    "-8 Southern"
    "-8 Waikato"
    "-8 Waikato"
    "-8 Waikato"
    "-8 Waitemata"
    "-8 Waitemata"
    "-8 Waitemata"
    "-8 Waitemata"
    "-8 Waitemata"
    "-8 Waitemata"
    "-8 Waitemata"
    "-9 Auckland"
    "-9 Auckland"
    "-9 Auckland"
    "-9 Bay.of.Plenty"
    "-9 Bay.of.Plenty"
    "-9 Canterbury"
    "-9 Canterbury"
    "-9 Canterbury"
    "-9 Canterbury"
    "-9 Canterbury"
    "-9 Canterbury"
    "-9 Capital.and.Coast"
    "-9 Counties.Manukau"
    "-9 Counties.Manukau"
    "-9 Counties.Manukau"
    "-9 Counties.Manukau"
    "-9 Hawkes.Bay"
    "-9 MidCentral"
    "-9 MidCentral"
    "-9 Northland"
    "-9 Southern"
    "-9 Southern"
    "-9 Southern"
    "-9 Southern"
    "-9 Southern"
    "-9 Southern"
    "-9 Southern"
    "-9 Southern"
    "-9 Waikato"
    "-9 Waikato"
    "-10 Auckland"
    "-10 Auckland"
    "-10 Auckland"
    "-10 Auckland"
    "-10 Auckland"
    "-10 Canterbury"
    "-10 Canterbury"
    "-10 Counties.Manukau"
    "-10 Counties.Manukau"
    "-10 Counties.Manukau"
    "-10 Counties.Manukau"
    "-10 Counties.Manukau"
    "-10 Hawkes.Bay"
    "-10 Hawkes.Bay"
    "-10 Hawkes.Bay"
    "-10 MidCentral"
    "-10 MidCentral"
    "-10 Nelson.Marlborough"
    "-10 Northland"
    "-10 Northland"
    "-10 Southern"
    "-10 Southern"
    "-10 Southern"
    "-10 Southern"
    "-10 Southern"
    "-10 Southern"
    "-10 Southern"
    "-10 Southern"
    "-10 Southern"
    "-10 Southern"
    "-10 Southern"
    "-10 Southern"
    "-10 Southern"
    "-10 Southern"
    "-10 Waikato"
    "-10 Waikato"
    "-10 Waikato"
    "-10 Waitemata"
    "-10 Waitemata"
    "-10 Waitemata"
    "-10 Waitemata"
    "-10 Waitemata"
    "-11 Auckland"
    "-11 Auckland"
    "-11 Auckland"
    "-11 Auckland"
    "-11 Auckland"
    "-11 Auckland"
    "-11 Auckland"
    "-11 Canterbury"
    "-11 Canterbury"
    "-11 Canterbury"
    "-11 Canterbury"
    "-11 Canterbury"
    "-11 Canterbury"
    "-11 Canterbury"
    "-11 Canterbury"
    "-11 Canterbury"
    "-11 Canterbury"
    "-11 Canterbury"
    "-11 Canterbury"
    "-11 Canterbury"
    "-11 Capital.and.Coast"
    "-11 Capital.and.Coast"
    "-11 Capital.and.Coast"
    "-11 Counties.Manukau"
    "-11 Counties.Manukau"
    "-11 Lakes"
    "-11 MidCentral"
    "-11 MidCentral"
    "-11 Nelson.Marlborough"
    "-11 Nelson.Marlborough"
    "-11 Nelson.Marlborough"
    "-11 Northland"
    "-11 Southern"
    "-11 Southern"
    "-11 Southern"
    "-11 Southern"
    "-11 Southern"
    "-11 Southern"
    "-11 Southern"
    "-11 Southern"
    "-11 Southern"
    "-11 Waikato"
    "-11 Waikato"
    "-11 Waikato"
    "-11 Waikato"
    "-11 Waitemata"
    "-11 Waitemata"
    "-11 Waitemata"
    "-11 Waitemata"
    "-12 Auckland"
    "-12 Auckland"
    "-12 Auckland"
    "-12 Canterbury"
    "-12 Canterbury"
    "-12 Canterbury"
    "-12 Canterbury"
    "-12 Canterbury"
    "-12 Capital.and.Coast"
    "-12 Capital.and.Coast"
    "-12 Capital.and.Coast"
    "-12 Counties.Manukau"
    "-12 Counties.Manukau"
    "-12 Counties.Manukau"
    "-12 Counties.Manukau"
    "-12 Counties.Manukau"
    "-12 Hawkes.Bay"
    "-12 Hawkes.Bay"
    "-12 Hawkes.Bay"
    "-12 Hawkes.Bay"
    "-12 Hawkes.Bay"
    "-12 Nelson.Marlborough"
    "-12 Northland"
    "-12 Northland"
    "-12 Southern"
    "-12 Southern"
    "-12 Southern"
    "-12 Southern"
    "-12 Southern"
    "-12 Southern"
    "-12 Southern"
    "-12 Southern"
    "-12 Southern"
    "-12 Southern"
    "-12 Southern"
    "-12 Southern"
    "-12 Southern"
    "-12 Southern"
    "-12 Waikato"
    "-12 Waikato"
    "-12 Waikato"
    "-12 Waikato"
    "-12 Waikato"
    "-12 Waikato"
    "-12 Waikato"
    "-12 Waikato"
    "-12 Waikato"
    "-12 Waikato"
    "-12 Waikato"
    "-12 Waitemata"
    "-12 Waitemata"
    "-12 Waitemata"
    "-12 Waitemata"
    "-12 Waitemata"
    "-12 Waitemata"
    "-13 Auckland"
    "-13 Auckland"
    "-13 Auckland"
    "-13 Auckland"
    "-13 Auckland"
    "-13 Bay.of.Plenty"
    "-13 Bay.of.Plenty"
    "-13 Bay.of.Plenty"
    "-13 Canterbury"
    "-13 Canterbury"
    "-13 Canterbury"
    "-13 Capital.and.Coast"
    "-13 Capital.and.Coast"
    "-13 Counties.Manukau"
    "-13 Counties.Manukau"
    "-13 Counties.Manukau"
    "-13 Counties.Manukau"
    "-13 Hawkes.Bay"
    "-13 Hawkes.Bay"
    "-13 Hutt.Valley"
    "-13 Hutt.Valley"
    "-13 Hutt.Valley"
    "-13 Hutt.Valley"
    "-13 MidCentral"
    "-13 Nelson.Marlborough"
    "-13 Nelson.Marlborough"
    "-13 Nelson.Marlborough"
    "-13 Northland"
    "-13 Southern"
    "-13 Southern"
    "-13 Southern"
    "-13 Southern"
    "-13 Southern"
    "-13 Southern"
    "-13 Southern"
    "-13 Southern"
    "-13 Southern"
    "-13 Southern"
    "-13 Southern"
    "-13 Southern"
    "-13 Southern"
    "-13 Southern"
    "-13 Southern"
    "-13 Southern"
    "-13 Southern"
    "-13 Waikato"
    "-13 Waikato"
    "-13 Waikato"
    "-13 Waikato"
    "-13 Waikato"
    "-13 Waikato"
    "-13 Waikato"
    "-13 Waitemata"
    "-13 Waitemata"
    "-13 Waitemata"
    "-13 Waitemata"
    "-13 Whanganui"
    "-13 Whanganui"
    "-14 Auckland"
    "-14 Auckland"
    "-14 Auckland"
    "-14 Auckland"
    "-14 Auckland"
    "-14 Auckland"
    "-14 Auckland"
    "-14 Bay.of.Plenty"
    "-14 Bay.of.Plenty"
    "-14 Bay.of.Plenty"
    "-14 Bay.of.Plenty"
    "-14 Bay.of.Plenty"
    "-14 Capital.and.Coast"
    "-14 Counties.Manukau"
    "-14 Counties.Manukau"
    "-14 Counties.Manukau"
    "-14 Hawkes.Bay"
    "-14 Hawkes.Bay"
    "-14 Hawkes.Bay"
    "-14 Hawkes.Bay"
    "-14 MidCentral"
    "-14 MidCentral"
    "-14 Nelson.Marlborough"
    "-14 Nelson.Marlborough"
    "-14 Nelson.Marlborough"
    "-14 Northland"
    "-14 Southern"
    "-14 Southern"
    "-14 Southern"
    "-14 Southern"
    "-14 Southern"
    "-14 Southern"
    "-14 Southern"
    "-14 Southern"
    "-14 Southern"
    "-14 Southern"
    "-14 Southern"
    "-14 Southern"
    "-14 Southern"
    "-14 Southern"
    "-14 Southern"
    "-14 Taranaki"
    "-14 Waikato"
    "-14 Waikato"
    "-14 Waikato"
    "-14 Waikato"
    "-14 Waikato"
    "-14 Waikato"
    "-14 Waikato"
    "-14 Waikato"
    "-14 Waikato"
    "-14 Waitemata"
    "-14 Waitemata"
    "-14 Waitemata"
    "-14 Waitemata"
    "-14 Waitemata"
    "-14 Waitemata"
    "-14 Waitemata"
    "-14 West.Coast"
    "-14 Whanganui"
    "-14 Whanganui"
    "-15 Auckland"
    "-15 Auckland"
    "-15 Auckland"
    "-15 Auckland"
    "-15 Bay.of.Plenty"
    "-15 Bay.of.Plenty"
    "-15 Canterbury"
    "-15 Canterbury"
    "-15 Capital.and.Coast"
    "-15 Capital.and.Coast"
    "-15 Counties.Manukau"
    "-15 Counties.Manukau"
    "-15 Counties.Manukau"
    "-15 Counties.Manukau"
    "-15 Counties.Manukau"
    "-15 Hutt.Valley"
    "-15 MidCentral"
    "-15 MidCentral"
    "-15 MidCentral"
    "-15 MidCentral"
    "-15 Nelson.Marlborough"
    "-15 Nelson.Marlborough"
    "-15 Northland"
    "-15 Southern"
    "-15 Southern"
    "-15 Southern"
    "-15 Southern"
    "-15 Southern"
    "-15 Southern"
    "-15 Southern"
    "-15 Waikato"
    "-15 Waikato"
    "-15 Waikato"
    "-15 Waikato"
    "-15 Waikato"
    "-15 Waikato"
    "-15 Waikato"
    "-15 Waikato"
    "-15 Waikato"
    "-15 Waikato"
    "-15 Waikato"
    "-15 Waikato"
    "-15 Waikato"
    "-15 Waikato"
    "-15 Waikato"
    "-15 Waikato"
    "-15 Waikato"
    "-15 Waikato"
    "-15 Waikato"
    "-15 Waikato"
    "-15 Waikato"
    "-15 Waikato"
    "-15 Waikato"
    "-15 Wairarapa"
    "-15 Waitemata"
    "-15 Waitemata"
    "-15 Waitemata"
    "-15 Waitemata"
    "-15 Waitemata"
    "-15 Waitemata"
    "-15 Waitemata"
    "-15 Waitemata"
    "-15 Waitemata"
    "-16 Auckland"
    "-16 Auckland"
    "-16 Auckland"
    "-16 Auckland"
    "-16 Bay.of.Plenty"
    "-16 Bay.of.Plenty"
    "-16 Bay.of.Plenty"
    "-16 Bay.of.Plenty"
    "-16 Canterbury"
    "-16 Counties.Manukau"
    "-16 Counties.Manukau"
    "-16 Counties.Manukau"
    "-16 Counties.Manukau"
    "-16 Counties.Manukau"
    "-16 Counties.Manukau"
    "-16 Hawkes.Bay"
    "-16 Hawkes.Bay"
    "-16 Hawkes.Bay"
    "-16 Lakes"
    "-16 Lakes"
    "-16 Nelson.Marlborough"
    "-16 Northland"
    "-16 South.Canterbury"
    "-16 Southern"
    "-16 Southern"
    "-16 Southern"
    "-16 Southern"
    "-16 Southern"
    "-16 Southern"
    "-16 Southern"
    "-16 Southern"
    "-16 Southern"
    "-16 Taranaki"
    "-16 Waikato"
    "-16 Waikato"
    "-16 Waikato"
    "-16 Waikato"
    "-16 Waitemata"
    "-16 Waitemata"
    "-16 Waitemata"
    "-16 Waitemata"
    "-16 Waitemata"
    "-17 Auckland"
    "-17 Auckland"
    "-17 Auckland"
    "-17 Auckland"
    "-17 Auckland"
    "-17 Auckland"
    "-17 Auckland"
    "-17 Auckland"
    "-17 Auckland"
    "-17 Canterbury"
    "-17 Canterbury"
    "-17 Capital.and.Coast"
    "-17 Capital.and.Coast"
    "-17 Capital.and.Coast"
    "-17 Capital.and.Coast"
    "-17 Capital.and.Coast"
    "-17 Capital.and.Coast"
    "-17 Capital.and.Coast"
    "-17 Counties.Manukau"
    "-17 Counties.Manukau"
    "-17 Counties.Manukau"
    "-17 MidCentral"
    "-17 MidCentral"
    "-17 Nelson.Marlborough"
    "-17 Northland"
    "-17 Northland"
    "-17 South.Canterbury"
    "-17 Southern"
    "-17 Southern"
    "-17 Southern"
    "-17 Southern"
    "-17 Southern"
    "-17 Southern"
    "-17 Southern"
    "-17 Southern"
    "-17 Southern"
    "-17 Southern"
    "-17 Southern"
    "-17 Southern"
    "-17 Waikato"
    "-17 Waikato"
    "-17 Waikato"
    "-17 Waikato"
    "-17 Waikato"
    "-17 Waikato"
    "-17 Waikato"
    "-17 Waikato"
    "-17 Waikato"
    "-17 Waikato"
    "-17 Waikato"
    "-17 Waikato"
    "-17 Waikato"
    "-17 Waikato"
    "-17 Waikato"
    "-17 Waikato"
    "-17 Waitemata"
    "-17 Waitemata"
    "-17 Waitemata"
    "-17 Waitemata"
    "-17 Waitemata"
    "-17 Waitemata"
    "-17 Waitemata"
    "-17 Waitemata"
    "-17 Waitemata"
    "-17 Waitemata"
    "-17 Waitemata"
    "-17 Waitemata"
    "-17 Waitemata"
    "-18 Auckland"
    "-18 Auckland"
    "-18 Auckland"
    "-18 Auckland"
    "-18 Auckland"
    "-18 Auckland"
    "-18 Auckland"
    "-18 Auckland"
    "-18 Auckland"
    "-18 Auckland"
    "-18 Auckland"
    "-18 Auckland"
    "-18 Auckland"
    "-18 Auckland"
    "-18 Bay.of.Plenty"
    "-18 Bay.of.Plenty"
    "-18 Canterbury"
    "-18 Canterbury"
    "-18 Canterbury"
    "-18 Canterbury"
    "-18 Capital.and.Coast"
    "-18 Capital.and.Coast"
    "-18 Capital.and.Coast"
    "-18 Capital.and.Coast"
    "-18 Capital.and.Coast"
    "-18 Counties.Manukau"
    "-18 Counties.Manukau"
    "-18 Counties.Manukau"
    "-18 Counties.Manukau"
    "-18 Counties.Manukau"
    "-18 Hawkes.Bay"
    "-18 Lakes"
    "-18 Nelson.Marlborough"
    "-18 South.Canterbury"
    "-18 Southern"
    "-18 Southern"
    "-18 Southern"
    "-18 Southern"
    "-18 Southern"
    "-18 Southern"
    "-18 Southern"
    "-18 Southern"
    "-18 Southern"
    "-18 Southern"
    "-18 Southern"
    "-18 Southern"
    "-18 Southern"
    "-18 Waikato"
    "-18 Waikato"
    "-18 Waikato"
    "-18 Waitemata"
    "-18 Waitemata"
    "-18 Waitemata"
    "-18 Waitemata"
    "-18 Waitemata"
    "-18 Waitemata"
    "-18 Waitemata"
    "-18 Waitemata"
    "-18 Waitemata"
    "-18 Waitemata"
    "-18 Waitemata"
    "-18 Waitemata"
    "-18 Waitemata"
    "-19 Auckland"
    "-19 Auckland"
    "-19 Auckland"
    "-19 Auckland"
    "-19 Auckland"
    "-19 Auckland"
    "-19 Auckland"
    "-19 Bay.of.Plenty"
    "-19 Bay.of.Plenty"
    "-19 Bay.of.Plenty"
    "-19 Canterbury"
    "-19 Canterbury"
    "-19 Canterbury"
    "-19 Capital.and.Coast"
    "-19 Capital.and.Coast"
    "-19 Capital.and.Coast"
    "-19 Counties.Manukau"
    "-19 Counties.Manukau"
    "-19 Counties.Manukau"
    "-19 Counties.Manukau"
    "-19 Hawkes.Bay"
    "-19 Hutt.Valley"
    "-19 Hutt.Valley"
    "-19 Hutt.Valley"
    "-19 MidCentral"
    "-19 MidCentral"
    "-19 South.Canterbury"
    "-19 South.Canterbury"
    "-19 South.Canterbury"
    "-19 South.Canterbury"
    "-19 Southern"
    "-19 Southern"
    "-19 Southern"
    "-19 Southern"
    "-19 Southern"
    "-19 Southern"
    "-19 Southern"
    "-19 Southern"
    "-19 Southern"
    "-19 Southern"
    "-19 Southern"
    "-19 Southern"
    "-19 Southern"
    "-19 Southern"
    "-19 Southern"
    "-19 Southern"
    "-19 Southern"
    "-19 Southern"
    "-19 Southern"
    "-19 Tairawhiti"
    "-19 Taranaki"
    "-19 Taranaki"
    "-19 Taranaki"
    "-19 Taranaki"
    "-19 Waikato"
    "-19 Waikato"
    "-19 Waikato"
    "-19 Waikato"
    "-19 Waikato"
    "-19 Waikato"
    "-19 Waikato"
    "-19 Waikato"
    "-19 Waikato"
    "-19 Waikato"
    "-19 Waikato"
    "-19 Waitemata"
    "-19 Waitemata"
    "-19 Whanganui"
    "-19 Whanganui"
    "-19 Whanganui"
    "-20 Auckland"
    "-20 Auckland"
    "-20 Auckland"
    "-20 Auckland"
    "-20 Auckland"
    "-20 Auckland"
    "-20 Auckland"
    "-20 Auckland"
    "-20 Auckland"
    "-20 Auckland"
    "-20 Auckland"
    "-20 Bay.of.Plenty"
    "-20 Bay.of.Plenty"
    "-20 Canterbury"
    "-20 Canterbury"
    "-20 Canterbury"
    "-20 Canterbury"
    "-20 Canterbury"
    "-20 Canterbury"
    "-20 Canterbury"
    "-20 Canterbury"
    "-20 Canterbury"
    "-20 Capital.and.Coast"
    "-20 Capital.and.Coast"
    "-20 Counties.Manukau"
    "-20 Counties.Manukau"
    "-20 Counties.Manukau"
    "-20 Hawkes.Bay"
    "-20 Hawkes.Bay"
    "-20 Hawkes.Bay"
    "-20 Hawkes.Bay"
    "-20 Hutt.Valley"
    "-20 Hutt.Valley"
    "-20 Hutt.Valley"
    "-20 Hutt.Valley"
    "-20 Lakes"
    "-20 Northland"
    "-20 Southern"
    "-20 Southern"
    "-20 Southern"
    "-20 Southern"
    "-20 Southern"
    "-20 Southern"
    "-20 Southern"
    "-20 Tairawhiti"
    "-20 Taranaki"
    "-20 Taranaki"
    "-20 Waikato"
    "-20 Waikato"
    "-20 Waikato"
    "-20 Waikato"
    "-20 Waikato"
    "-20 Waikato"
    "-20 Waikato"
    "-20 Waikato"
    "-20 Waikato"
    "-20 Waikato"
    "-20 Waikato"
    "-20 Waitemata"
    "-20 Waitemata"
    "-20 Waitemata"
    "-20 Waitemata"
    "-20 Waitemata"
    "-20 Waitemata"
    "-20 Waitemata"
    "-20 Waitemata"
    "-20 West.Coast"
    "-20 West.Coast"
    "-21 Auckland"
    "-21 Auckland"
    "-21 Auckland"
    "-21 Auckland"
    "-21 Auckland"
    "-21 Auckland"
    "-21 Auckland"
    "-21 Auckland"
    "-21 Auckland"
    "-21 Auckland"
    "-21 Auckland"
    "-21 Bay.of.Plenty"
    "-21 Bay.of.Plenty"
    "-21 Canterbury"
    "-21 Canterbury"
    "-21 Canterbury"
    "-21 Canterbury"
    "-21 Capital.and.Coast"
    "-21 Capital.and.Coast"
    "-21 Capital.and.Coast"
    "-21 Capital.and.Coast"
    "-21 Counties.Manukau"
    "-21 Counties.Manukau"
    "-21 Counties.Manukau"
    "-21 Counties.Manukau"
    "-21 Counties.Manukau"
    "-21 Counties.Manukau"
    "-21 Hawkes.Bay"
    "-21 Hawkes.Bay"
    "-21 Hawkes.Bay"
    "-21 Hutt.Valley"
    "-21 Lakes"
    "-21 Lakes"
    "-21 Lakes"
    "-21 MidCentral"
    "-21 Nelson.Marlborough"
    "-21 Nelson.Marlborough"
    "-21 Nelson.Marlborough"
    "-21 Northland"
    "-21 Northland"
    "-21 Southern"
    "-21 Southern"
    "-21 Southern"
    "-21 Southern"
    "-21 Southern"
    "-21 Southern"
    "-21 Southern"
    "-21 Southern"
    "-21 Southern"
    "-21 Taranaki"
    "-21 Waikato"
    "-21 Waikato"
    "-21 Waikato"
    "-21 Waikato"
    "-21 Waikato"
    "-21 Waikato"
    "-21 Waikato"
    "-21 Waikato"
    "-21 Waikato"
    "-21 Waikato"
    "-21 Waikato"
    "-21 Waikato"
    "-21 Waikato"
    "-21 Waikato"
    "-21 Waikato"
    "-21 Waikato"
    "-21 Waikato"
    "-21 Waikato"
    "-21 Waikato"
    "-21 Waikato"
    "-21 Waitemata"
    "-21 Waitemata"
    "-21 Waitemata"
    "-21 Waitemata"
    "-21 Waitemata"
    "-21 Waitemata"
    "-21 Waitemata"
    "-21 Waitemata"
    "-21 Waitemata"
    "-21 Waitemata"
    "-22 Auckland"
    "-22 Auckland"
    "-22 Auckland"
    "-22 Auckland"
    "-22 Auckland"
    "-22 Auckland"
    "-22 Bay.of.Plenty"
    "-22 Canterbury"
    "-22 Canterbury"
    "-22 Canterbury"
    "-22 Canterbury"
    "-22 Canterbury"
    "-22 Canterbury"
    "-22 Canterbury"
    "-22 Capital.and.Coast"
    "-22 Capital.and.Coast"
    "-22 Capital.and.Coast"
    "-22 Capital.and.Coast"
    "-22 Capital.and.Coast"
    "-22 Capital.and.Coast"
    "-22 Capital.and.Coast"
    "-22 Capital.and.Coast"
    "-22 Capital.and.Coast"
    "-22 Capital.and.Coast"
    "-22 Counties.Manukau"
    "-22 Counties.Manukau"
    "-22 Hawkes.Bay"
    "-22 Hawkes.Bay"
    "-22 Lakes"
    "-22 MidCentral"
    "-22 Nelson.Marlborough"
    "-22 Nelson.Marlborough"
    "-22 Northland"
    "-22 Northland"
    "-22 Southern"
    "-22 Southern"
    "-22 Southern"
    "-22 Southern"
    "-22 Southern"
    "-22 Southern"
    "-22 Southern"
    "-22 Waikato"
    "-22 Waikato"
    "-22 Wairarapa"
    "-22 Wairarapa"
    "-22 Waitemata"
    "-22 Waitemata"
    "-22 Waitemata"
    "-22 Waitemata"
    "-22 Waitemata"
    "-22 Waitemata"
    "-22 Waitemata"
    "-22 Waitemata"
    "-22 Waitemata"
    "-23 Auckland"
    "-23 Auckland"
    "-23 Auckland"
    "-23 Auckland"
    "-23 Auckland"
    "-23 Auckland"
    "-23 Auckland"
    "-23 Canterbury"
    "-23 Canterbury"
    "-23 Canterbury"
    "-23 Canterbury"
    "-23 Capital.and.Coast"
    "-23 Capital.and.Coast"
    "-23 Capital.and.Coast"
    "-23 Capital.and.Coast"
    "-23 Capital.and.Coast"
    "-23 Capital.and.Coast"
    "-23 Capital.and.Coast"
    "-23 Capital.and.Coast"
    "-23 Counties.Manukau"
    "-23 Hutt.Valley"
    "-23 Lakes"
    "-23 MidCentral"
    "-23 Nelson.Marlborough"
    "-23 Nelson.Marlborough"
    "-23 Nelson.Marlborough"
    "-23 South.Canterbury"
    "-23 South.Canterbury"
    "-23 Southern"
    "-23 Southern"
    "-23 Southern"
    "-23 Southern"
    "-23 Waikato"
    "-23 Waikato"
    "-23 Waikato"
    "-23 Waikato"
    "-23 Waitemata"
    "-23 Waitemata"
    "-24 Auckland"
    "-24 Auckland"
    "-24 Auckland"
    "-24 Auckland"
    "-24 Auckland"
    "-24 Auckland"
    "-24 Auckland"
    "-24 Canterbury"
    "-24 Canterbury"
    "-24 Canterbury"
    "-24 Capital.and.Coast"
    "-24 Capital.and.Coast"
    "-24 Capital.and.Coast"
    "-24 Capital.and.Coast"
    "-24 Capital.and.Coast"
    "-24 Counties.Manukau"
    "-24 Counties.Manukau"
    "-24 Counties.Manukau"
    "-24 Counties.Manukau"
    "-24 Counties.Manukau"
    "-24 Counties.Manukau"
    "-24 Counties.Manukau"
    "-24 Counties.Manukau"
    "-24 Hutt.Valley"
    "-24 Nelson.Marlborough"
    "-24 Nelson.Marlborough"
    "-24 Nelson.Marlborough"
    "-24 Nelson.Marlborough"
    "-24 Nelson.Marlborough"
    "-24 Northland"
    "-24 Southern"
    "-24 Southern"
    "-24 Waikato"
    "-24 Waikato"
    "-24 Waitemata"
    "-24 Waitemata"
    "-24 Waitemata"
    "-24 Waitemata"
    "-24 Waitemata"
    "-24 West.Coast"
    "-25 Auckland"
    "-25 Auckland"
    "-25 Auckland"
    "-25 Auckland"
    "-25 Bay.of.Plenty"
    "-25 Canterbury"
    "-25 Canterbury"
    "-25 Canterbury"
    "-25 Capital.and.Coast"
    "-25 Counties.Manukau"
    "-25 Counties.Manukau"
    "-25 Counties.Manukau"
    "-25 Counties.Manukau"
    "-25 MidCentral"
    "-25 Nelson.Marlborough"
    "-25 Nelson.Marlborough"
    "-25 Taranaki"
    "-25 Taranaki"
    "-25 Waikato"
    "-25 Waikato"
    "-25 Waikato"
    "-25 Wairarapa"
    "-25 Wairarapa"
    "-25 Waitemata"
    "-26 Auckland"
    "-26 Canterbury"
    "-26 Canterbury"
    "-26 Lakes"
    "-26 Nelson.Marlborough"
    "-26 Nelson.Marlborough"
    "-26 Southern"
    "-26 Southern"
    "-26 Waikato"
    "-26 Waitemata"
    "-26 Waitemata"
    "-27 Auckland"
    "-27 Auckland"
    "-27 Capital.and.Coast"
    "-27 Capital.and.Coast"
    "-27 Counties.Manukau"
    "-27 Hawkes.Bay"
    "-27 Hutt.Valley"
    "-27 Hutt.Valley"
    "-27 MidCentral"
    "-27 MidCentral"
    "-27 MidCentral"
    "-27 Southern"
    "-27 Waikato"
    "-27 Waikato"
    "-27 Wairarapa"
    "-27 Waitemata"
    "-28 Auckland"
    "-28 Canterbury"
    "-28 Capital.and.Coast"
    "-28 Southern"
    "-28 Southern"
    "-29 Auckland"
    "-29 Auckland"
    "-29 Canterbury"
    "-29 Capital.and.Coast"
    "-29 Counties.Manukau"
    "-29 Lakes"
    "-29 Northland"
    "-29 Southern"
    "-29 Southern"
    "-29 Waitemata"
    "-29 Waitemata"
    "-30 Auckland"
    "-30 Capital.and.Coast"
    "-30 Capital.and.Coast"
    "-30 Taranaki"
    "-30 Taranaki"
    "-30 Waikato"
    "-30 Waikato"
    "-31 Waitemata"
    "-32 Capital.and.Coast"
    "-33 Southern"
    "-34 Counties.Manukau"
    "-40 Counties.Manukau"
    "-42 Counties.Manukau"
    "-42 Waitemata"
    "-44 Waitemata"
    "-49 Auckland"
    "-1 Auckland"
    "-1 Bay.of.Plenty"
    "-1 Bay.of.Plenty"
    "-1 Bay.of.Plenty"
    "-1 Counties.Manukau"
    "-1 Southern"
    "-1 Southern"
    "-1 Waikato"
    "-2 Auckland"
    "-2 Southern"
    "-2 Waitemata"
    "-2 Waitemata"
    "-2 Waitemata"
    "-2 Waitemata"
    "-2 Waitemata"
    "-3 Auckland"
    "-3 Auckland"
    "-3 Auckland"
    "-3 Canterbury"
    "-3 Capital.and.Coast"
    "-3 Counties.Manukau"
    "-3 Southern"
    "-3 Southern"
    "-3 Southern"
    "-3 Southern"
    "-4 Auckland"
    "-4 Auckland"
    "-4 Canterbury"
    "-4 Waikato"
    "-4 Waikato"
    "-4 Waikato"
    "-4 Waikato"
    "-4 Waikato"
    "-4 Waikato"
    "-4 Waitemata"
    "-4 Waitemata"
    "-4 Waitemata"
    "-5 Auckland"
    "-5 Auckland"
    "-5 Auckland"
    "-5 Auckland"
    "-5 Bay.of.Plenty"
    "-5 Canterbury"
    "-5 Canterbury"
    "-5 Hawkes.Bay"
    "-5 Lakes"
    "-5 Southern"
    "-5 Southern"
    "-5 Southern"
    "-5 Southern"
    "-5 Southern"
    "-5 Waikato"
    "-5 Waikato"
    "-5 Waitemata"
    "-5 Waitemata"
    "-5 Waitemata"
    "-5 Waitemata"
    "-6 Auckland"
    "-6 Auckland"
    "-6 Auckland"
    "-6 Bay.of.Plenty"
    "-6 Bay.of.Plenty"
    "-6 Canterbury"
    "-6 Canterbury"
    "-6 Canterbury"
    "-6 Canterbury"
    "-6 Canterbury"
    "-6 Hawkes.Bay"
    "-6 Northland"
    "-6 Southern"
    "-6 Southern"
    "-7 Auckland"
    "-7 Canterbury"
    "-7 Canterbury"
    "-7 Canterbury"
    "-7 Southern"
    "-7 Waikato"
    "-7 Waikato"
    "-7 Waikato"
    "-7 Waitemata"
    "-8 Auckland"
    "-8 Auckland"
    "-8 Canterbury"
    "-8 Canterbury"
    "-8 Canterbury"
    "-8 Counties.Manukau"
    "-8 Nelson.Marlborough"
    "-8 Nelson.Marlborough"
    "-8 Southern"
    "-8 Southern"
    "-8 Waikato"
    "-8 Waikato"
    "-8 Waitemata"
    "-8 Waitemata"
    "-9 Auckland"
    "-9 Auckland"
    "-9 Auckland"
    "-9 Bay.of.Plenty"
    "-9 Bay.of.Plenty"
    "-9 Canterbury"
    "-9 Canterbury"
    "-9 Canterbury"
    "-9 Canterbury"
    "-9 Canterbury"
    "-9 Canterbury"
    "-9 Capital.and.Coast"
    "-9 Counties.Manukau"
    "-9 Counties.Manukau"
    "-9 Counties.Manukau"
    "-9 Nelson.Marlborough"
    "-9 Nelson.Marlborough"
    "-9 Waikato"
    "-9 Waikato"
    "-9 Waikato"
    "-9 Waikato"
    "-9 Waikato"
    "-9 Waitemata"
    "-10 Auckland"
    "-10 Auckland"
    "-10 Auckland"
    "-10 Auckland"
    "-10 Auckland"
    "-10 Bay.of.Plenty"
    "-10 Canterbury"
    "-10 Canterbury"
    "-10 Canterbury"
    "-10 Counties.Manukau"
    "-10 Counties.Manukau"
    "-10 Counties.Manukau"
    "-10 Hawkes.Bay"
    "-10 Hutt.Valley"
    "-10 MidCentral"
    "-10 Nelson.Marlborough"
    "-10 Nelson.Marlborough"
    "-10 Nelson.Marlborough"
    "-10 Nelson.Marlborough"
    "-10 Southern"
    "-10 Southern"
    "-10 Southern"
    "-10 Waikato"
    "-10 Waitemata"
    "-10 Waitemata"
    "-10 Waitemata"
    "-10 Waitemata"
    "-10 West.Coast"
    "-11 Auckland"
    "-11 Auckland"
    "-11 Auckland"
    "-11 Auckland"
    "-11 Canterbury"
    "-11 Canterbury"
    "-11 Canterbury"
    "-11 Canterbury"
    "-11 Canterbury"
    "-11 Canterbury"
    "-11 Capital.and.Coast"
    "-11 Capital.and.Coast"
    "-11 Capital.and.Coast"
    "-11 Counties.Manukau"
    "-11 MidCentral"
    "-11 Nelson.Marlborough"
    "-11 Nelson.Marlborough"
    "-11 Southern"
    "-11 Waikato"
    "-11 Waikato"
    "-11 Waitemata"
    "-11 Waitemata"
    "-11 Waitemata"
    "-11 Waitemata"
    "-11 Waitemata"
    "-11 Waitemata"
    "-11 Waitemata"
    "-11 Waitemata"
    "-11 Waitemata"
    "-11 Waitemata"
    "-12 Auckland"
    "-12 Auckland"
    "-12 Auckland"
    "-12 Bay.of.Plenty"
    "-12 Canterbury"
    "-12 Canterbury"
    "-12 Canterbury"
    "-12 Canterbury"
    "-12 Canterbury"
    "-12 Capital.and.Coast"
    "-12 Capital.and.Coast"
    "-12 Capital.and.Coast"
    "-12 Hutt.Valley"
    "-12 Southern"
    "-12 Waikato"
    "-12 Waikato"
    "-12 Waitemata"
    "-12 Waitemata"
    "-13 Bay.of.Plenty"
    "-13 Canterbury"
    "-13 Canterbury"
    "-13 Canterbury"
    "-13 Canterbury"
    "-13 Canterbury"
    "-13 Canterbury"
    "-13 Capital.and.Coast"
    "-13 Capital.and.Coast"
    "-13 Counties.Manukau"
    "-13 Counties.Manukau"
    "-13 Hawkes.Bay"
    "-13 Hawkes.Bay"
    "-13 Hawkes.Bay"
    "-13 Southern"
    "-13 Waikato"
    "-13 Wairarapa"
    "-13 Wairarapa"
    "-13 Waitemata"
    "-13 Waitemata"
    "-13 Waitemata"
    "-13 Waitemata"
    "-13 Waitemata"
    "-13 Waitemata"
    "-13 Waitemata"
    "-14 Auckland"
    "-14 Canterbury"
    "-14 Canterbury"
    "-14 Capital.and.Coast"
    "-14 Capital.and.Coast"
    "-14 Nelson.Marlborough"
    "-14 Southern"
    "-14 Waikato"
    "-14 Waitemata"
    "-14 Waitemata"
    "-14 Waitemata"
    "-14 Waitemata"
    "-14 Waitemata"
    "-15 Auckland"
    "-15 Auckland"
    "-15 Auckland"
    "-15 Auckland"
    "-15 Bay.of.Plenty"
    "-15 Canterbury"
    "-15 Counties.Manukau"
    "-15 Counties.Manukau"
    "-15 Southern"
    "-15 Southern"
    "-15 Waikato"
    "-15 Waitemata"
    "-15 Waitemata"
    "-15 Waitemata"
    "-16 Auckland"
    "-16 Capital.and.Coast"
    "-16 Counties.Manukau"
    "-16 Counties.Manukau"
    "-16 Nelson.Marlborough"
    "-16 Waikato"
    "-16 Waitemata"
    "-16 Waitemata"
    "-16 Waitemata"
    "-16 Waitemata"
    "-17 Auckland"
    "-17 Canterbury"
    "-17 Capital.and.Coast"
    "-17 MidCentral"
    "-17 South.Canterbury"
    "-17 Waitemata"
    "-17 Waitemata"
    "-18 Auckland"
    "-18 Auckland"
    "-18 Canterbury"
    "-18 Canterbury"
    "-18 Canterbury"
    "-18 Canterbury"
    "-18 Counties.Manukau"
    "-19 Auckland"
    "-19 Canterbury"
    "-19 Canterbury"
    "-19 Canterbury"
    "-19 Waitemata"
    "-20 Auckland"
    "-20 Bay.of.Plenty"
    "-20 Canterbury"
    "-20 Canterbury"
    "-20 Canterbury"
    "-20 Capital.and.Coast"
    "-20 Nelson.Marlborough"
    "-20 Northland"
    "-20 Taranaki"
    "-20 Waitemata"
    "-20 Waitemata"
    "-20 Waitemata"
    "-21 Canterbury"
    "-21 MidCentral"
    "-21 Waikato"
    "-21 Waitemata"
    "-22 Bay.of.Plenty"
    "-22 Bay.of.Plenty"
    "-22 Nelson.Marlborough"
    "-22 Nelson.Marlborough"
    "-22 Waitemata"
    "-23 Canterbury"
    "-23 Capital.and.Coast"
    "-23 Hutt.Valley"
    "-24 Auckland"
    "-24 Canterbury"
    "-24 Capital.and.Coast"
    "-24 Capital.and.Coast"
    "-24 Waikato"
    "-27 Auckland"
    "-27 Auckland"
    "-29 Southern"
    "-41 Counties.Manukau"
    "-41 Waitemata" ] "\n"
end
@#$#@#$#@
GRAPHICS-WINDOW
533
12
1034
746
-1
-1
29.0
1
12
1
1
1
0
0
0
1
-8
8
-12
12
1
1
1
days
100.0

BUTTON
426
68
517
103
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
425
107
518
142
go-one-day
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

BUTTON
422
255
518
290
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
149
183
182
num-locales
num-locales
20
200
20.0
10
1
NIL
HORIZONTAL

SLIDER
1409
479
1549
512
exp-to-presymp
exp-to-presymp
0
1
0.2
0.01
1
NIL
HORIZONTAL

SLIDER
1410
519
1548
552
presymp-to-inf
presymp-to-inf
0
1
0.2
0.01
1
NIL
HORIZONTAL

SLIDER
1409
559
1547
592
inf-to-rec
inf-to-rec
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
1407
599
1545
632
rel-inf-presymp
rel-inf-presymp
0
1
0.15
0.01
1
NIL
HORIZONTAL

MONITOR
1425
404
1532
449
mean-trans-coeff
mean-trans-coeff
5
1
11

SLIDER
1405
740
1545
773
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
1405
779
1547
812
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
1406
660
1544
693
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
1406
700
1545
733
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
1404
819
1547
852
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
262
69
406
102
initial-infected
initial-infected
0
10000
600.0
10
1
NIL
HORIZONTAL

PLOT
1049
12
1394
394
Current totals
days
log (people + 1)
0.0
6.0
0.0
3.0
true
true
"" ""
PENS
"exposed" 1.0 0 -8431303 true "" "plot log (sum [exposed] of locales + 1) 10"
"presymp" 1.0 0 -955883 true "" "plot log (sum [presymptomatic] of locales + 1) 10"
"infected" 1.0 0 -2674135 true "" "plot log (sum [infected] of locales + 1) 10"
"recovered" 1.0 0 -13840069 true "" "plot log (sum [recovered] of locales + 1) 10"
"dead" 1.0 0 -16777216 true "" "plot log (total-dead + 1) 10"

SLIDER
262
249
404
282
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
264
211
403
244
use-seed?
use-seed?
1
1
-1000

SLIDER
9
112
182
145
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
12
63
209
107
Population\nOnly applicable when not initialising from data
12
0.0
1

TEXTBOX
1406
644
1490
662
Mortality
12
0.0
1

TEXTBOX
8
292
92
310
Connectivity
12
0.0
1

TEXTBOX
9
474
141
493
Control and testing
12
0.0
1

SLIDER
9
186
181
219
pop-sd-multiplier
pop-sd-multiplier
0.01
1.2
0.75
0.01
1
NIL
HORIZONTAL

MONITOR
98
230
178
275
max-pop
max [pop-0] of locales
0
1
11

MONITOR
12
232
93
277
min-pop
min [pop-0] of locales
0
1
11

INPUTBOX
5
492
291
602
alert-levels-control
pessimistic [1 0.8 0.6 0.36]\nrealistic [1 0.72 0.52 0.32]\noptimistic [1 0.64 0.44 0.28]\nother [1 0.7 0.4 0.16]
1
1
String

MONITOR
1424
353
1531
398
eff-mean-R0
mean-R0
3
1
11

SLIDER
349
389
522
422
initial-alert-level
initial-alert-level
1
4
2.0
1
1
NIL
HORIZONTAL

MONITOR
1400
211
1473
256
dead
total-dead
0
1
11

PLOT
1048
402
1393
787
Daily counts
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
6
750
166
783
pop-test-rate
pop-test-rate
0
0.002
0.001
0.0001
1
NIL
HORIZONTAL

INPUTBOX
7
352
163
412
alert-levels-flow
[0.2 0.1 0.05 0.025]
1
0
String

CHOOSER
192
376
345
421
alert-policy
alert-policy
"static" "local" "global-mean" "global-max" "local-random"
0

MONITOR
1398
12
1524
57
total-infected
total-infected
0
1
11

MONITOR
1400
162
1473
207
recovered
total-recovered
0
1
11

INPUTBOX
314
534
521
594
alert-level-triggers
[0.0001 0.00025 0.0005 1]
1
0
String

SLIDER
318
637
520
670
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
530
760
701
793
log-all-locales?
log-all-locales?
1
1
-1000

SLIDER
320
598
522
631
start-lifting-quarantine
start-lifting-quarantine
0
56
7.0
7
1
days
HORIZONTAL

TEXTBOX
351
370
430
388
Alert levels
12
0.0
1

MONITOR
269
680
327
725
pop-lev-1
sum [pop-0] of locales with [alert-level = 1]
0
1
11

MONITOR
332
680
390
725
pop-lev-2
sum [pop-0] of locales with [alert-level = 2]
0
1
11

MONITOR
395
680
454
725
pop-lev-3
sum [pop-0] of locales with [alert-level = 3]
0
1
11

MONITOR
459
679
518
724
pop-lev-4
sum [pop-0] of locales with [alert-level = 4]
0
1
11

MONITOR
269
731
327
776
n-lev-1
count locales with [alert-level = 1]
0
1
11

MONITOR
332
731
390
776
n-lev-2
count locales with [alert-level = 2]
0
1
11

MONITOR
395
731
454
776
n-lev-3
count locales with [alert-level = 3]
0
1
11

MONITOR
459
730
518
775
n-lev-4
count locales with [alert-level = 4]
0
1
11

INPUTBOX
529
796
714
856
log-folder
foo
1
0
String

MONITOR
415
790
509
835
alert-activity
alert-level-changes / count locales
4
1
11

BUTTON
744
755
900
790
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
906
755
1032
789
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

TEXTBOX
1415
294
1515
312
Pandemic
12
0.0
1

TEXTBOX
11
842
481
885
Model parameters not tuned to any specific location.\nExercise caution in using model to inform decision making.
16
15.0
1

CHOOSER
9
10
197
55
setup-method
setup-method
"NZ DHBs random cases" "NZ TAs random cases" "NZ DHBs from Apr 15 MoH data" "Random landscape" "Costa Rica"
1

SLIDER
6
312
215
345
max-connection-distance
max-connection-distance
150
1200
600.0
25
1
km
HORIZONTAL

BUTTON
424
144
518
178
go-one-week
repeat 7 [go]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

BUTTON
423
181
518
215
go-4-weeks
repeat 28 [go]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

BUTTON
424
217
519
251
go-13-weeks
repeat 91 [go]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

BUTTON
422
294
518
328
reset
clear-all
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
222
313
406
346
new-exposures-arriving
new-exposures-arriving
0
50
0.0
0.1
1
NIL
HORIZONTAL

TEXTBOX
277
179
406
209
Use seed value to repeat a run exactly
12
0.0
1

TEXTBOX
223
298
388
317
Border control not total
12
0.0
1

CHOOSER
6
609
149
654
control-scenario
control-scenario
"pessimistic" "realistic" "optimistic" "other"
3

TEXTBOX
335
792
408
837
Alert level changes per locale
12
0.0
1

TEXTBOX
309
434
529
480
NOTE: with alert-policy 'static'\ninteractively change global level\nusing the initial-alert-level control
12
0.0
1

SLIDER
1412
312
1537
345
R0
R0
1
10
2.5
0.01
1
NIL
HORIZONTAL

TEXTBOX
1410
457
1547
476
Disease progression
12
0.0
1

OUTPUT
206
7
528
57
12

TEXTBOX
206
109
408
128
If not initialising from case data
12
0.0
1

MONITOR
171
739
253
784
daily-tests
sum [item 0 recent-tests] of locales
0
1
11

SWITCH
7
417
163
450
gravity-weight?
gravity-weight?
0
1
-1000

SWITCH
911
802
1033
835
debug?
debug?
1
1
-1000

SLIDER
8
710
167
743
inf-test-rate
inf-test-rate
0
1
0.9
0.01
1
NIL
HORIZONTAL

SLIDER
6
789
200
822
false-negative-rate
false-negative-rate
0
0.2
0.05
0.01
1
NIL
HORIZONTAL

SWITCH
228
131
406
164
initialise-by-burn-in?
initialise-by-burn-in?
0
1
-1000

MONITOR
1399
62
1525
107
NIL
total-presymptomatic
0
1
11

MONITOR
1399
112
1524
157
NIL
total-exposed
0
1
11

SLIDER
332
489
522
522
response-time
response-time
0
28
7.0
1
1
days
HORIZONTAL

INPUTBOX
215
590
320
665
script
0 4\n35 3\n56 2
1
1
String

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

circle 3
false
0
Circle -7500403 false true -1 -1 301
Circle -7500403 false true 14 14 272
Circle -7500403 false true 7 7 286

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
