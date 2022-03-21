;v2.2.0 has minor changes for better functionality.
extensions [ gis ]
globals
[
  dxp                            ;x coordinate of the patch where a yearling deer belongs before dispersal (natal range)
  dxn                            ;x coordinate of the patch where a yearling deer reaches after dispersing from its natal range
  dyp                            ;y coordinate of the patch where a yearling deer belongs before dispersal (natal range)
  dyn                            ;y coordinate of the patch where a yearling deer reaches after dispersing from its natal range
  dd                             ;predicted dispersal distance from log-normal pdf
  ndd                            ;counter-number of deer dispersing out of the model landscape
  vals1                          ;list to store reporters for output file
  vals2                          ;list to store reporters for output file
  vals3                          ;list to store reporters for output file
  vals                           ;list to store reporters for output file
  tmfh                           ;counter-male fawn deer harvested
  tmyh                           ;counter-male yearling deer harvested
  tamh                           ;counter-male adult deer harvested
  tffh                           ;counter-female fawn deer harvested
  tfyh                           ;counter-female yearling deer harvested
  tafh                           ;counter-female adult deer harvested
  tgroid                         ;stores groid (group-id) of an individual during implementation of certain submodels
  ttgroid                        ;stores groid (group-id) of an individual during implementation of certain submodels
  n_leaders_lost                 ;counter-doe social group leaders losing their leadership status
  twho                           ;stores 'who' of an individual during implementation of certain submodels
  counter1                       ;counter for fawns less than 2 month old while implementing hunting or non-hunting mortality for a female deer
  tgr                            ;counter for group members during group formation as well as fission
  oldm                           ;proportion of old males (above 229 when d = 1) in the adult male deer population
  oldf                           ;proportion of old females (above 229 when d = 1) in the adult female deer population
  tmgroid                        ;stores mgroid (male group id) during execution of certain submodels
  nt
  ntplus1
  lambda
  iteration
  ]
patches-own
[
  forest-percent                ;percent forest cover
  dh                            ;'deer habitat' ≥ 1 if a patch or a group of contiguous patches qualifies as deer habitat, otherwise < 1
  do                            ;'deer occupancy' 1 if deer occur on a patch, 0 otherwise
  border                        ;identifying border patches
  dfp                           ;mean forest-percent for a patch and its immediate neighbors
  ]

breed [ deers deer ]
deers-own
[
  sex                          ;gender of the individual
  aim                          ;age in months
  momid                        ;mother's who number
  gl                           ;1 if doe social group leader; 0 for group members and others
  groid                        ;group id, who of the doe social group leader; -1 if not a member of any doe social group
  gr                           ;takes the value equal to group size (not including the leader) for doe social group leaders; -1 for non-leader group member, and -2 for deer not in any group
  ml                           ;male bachelor group leader ml = 1 otherwise 0
  mgroid                       ;male bachelor group identifier: -2 at birth, -1 after dispersal and = who of leader when join group as adults
  ]
;-----------------------------------------------------------------------------------;
to setup
  clear-all
  setup-landscape
  if recommended_parameter_values = TRUE [
    ;demographic parameteres for Missouri deer populations (extrapolated from surveys/expert opinions)
    set sexratio 120
    set adultprop 0.40
    set yearlingprop 0.25
    set mf6nhm 0.055
    set ff6nhm 0.055
    set mf12nhm 0.05
    set ff12nhm 0.05
    set mynhm 0.01
    set fynhm 0.00
    set manhm 0.01
    set fanhm 0.02
    set mf6hm 0
    set ff6hm 0
    set mf12hm 0.05
    set ff12hm 0.02
    set myhm 0.25
    set fyhm 0.15
    set mahm 0.40
    set fahm 0.20 ]
  let ini-tot (mfd + myd + mad + ffd + fyd + fad)
  let ivals (list (mfd) (myd) (mad) (ffd) (fyd) (fad) (ini-tot))
  if file-exists? (word "../results/deerpopdy" region ".csv") = FALSE [
    file-open (word "../results/deerpopdy" region ".csv")
    file-print "year,posth_mf,posth_my,posth_ma,posth_ff,posth_fy,posth_fa,posth_total,preh_mf,preh_my,preh_ma,preh_ff,preh_fy,preh_fa,preh_total,lambda,mf_harvest,my_harvest,ma_harvest,ff_harvest,fy_harvest,fa_harvest,total_harvest,iteration"
    file-type first ivals
    foreach but-first ivals [ [ ?1 ] ->
      file-type "," file-type ?1
      ]
    file-print""
    file-close
    ]
  ;for sensitivity analysis
  if file-exists? (word "../results/sa" region ".csv") = FALSE [
    file-open (word "../results/sa" region ".csv")
    file-print "year,pre-harvest_total,adult_prop,yearling_prop,fawn_prop,f:m_ratio,lambda,iteration"
    file-close
    ]
  reset-ticks
end
;-----------------------------------------------------------------------------------;
to setup-landscape      ;setup model landscape using GIS data for the selected region
  let forest-map []
  if region = "StJoseph" [
    resize-world 0 25 0 23
    set forest-map gis:load-dataset "../data/stjoseph.asc"
    if recommended_parameter_values = TRUE [ set post_harvest_density 25 ]
    ]
  if region = "Elkhart" [
    resize-world 0 22 0 23
    set forest-map gis:load-dataset "../data/elkhart.asc"
    if recommended_parameter_values = TRUE [ set post_harvest_density 25 ]
    ]
  if region = "LaGrange" [
    resize-world 0 25 0 17
    set forest-map gis:load-dataset "../data/lagrange.asc"
    if recommended_parameter_values = TRUE [ set post_harvest_density 25 ]
    ]
  if region = "LaPorte" [
    resize-world 0 23 0 37
    set forest-map gis:load-dataset "../data/laporte.asc"
    if recommended_parameter_values = TRUE [ set post_harvest_density 25 ]
    ]
  if region = "Kankakee" [
    resize-world 0 122 0 108
    set forest-map gis:load-dataset "../data/kankakee_20220321.asc"
    if recommended_parameter_values = TRUE [ set post_harvest_density 25 ]
    ]

  gis:set-world-envelope gis:envelope-of forest-map
  gis:apply-raster forest-map forest-percent

  ask patches [ set border 0 ]
  let brown-patches patches with [ forest-percent = 0 ]
  ask brown-patches [
    set pcolor brown
    set dh -1
    ]
; patch procedure: identify border patches
  let pborder1 patches with-max [ pxcor ]
  ask pborder1 [ set border 1 ]
  let pborder2 patches with-min [ pxcor ]
  ask pborder2 [ set border 2 ]
  let pborder3 patches with-max [ pycor ]
  ask pborder3 [
    ifelse border > 0
    [ ifelse border > 1
      [ set border 8 ]                  ;upper left corner patch
      [ set border 5 ]                  ;upper right corner patch
    ]
    [ set border 3 ]
    ]
  let pborder4 patches with-min [ pycor ]
  ask pborder4 [
    ifelse border > 0
    [ifelse border > 1
      [ set border 7 ]                   ;lower left corner patch
      [ set border 6 ]                   ;lower right corner patch
      ]
    [ set border 4 ]
    ]
;patch procedure: identify 'deer occupancy' patches and 'deer habitat' patches
  let hrng 0
  let optfc patches with [ forest-percent >= min-forestcover-percent ]
  ask optfc [
    ifelse forest-percent <= max-forestcover-percent
    [ set do 1
      set hrng hrng + 1
      set dh hrng
      ]
    [ set do 2 ]
    ]
  let belowoptfc patches with [ forest-percent < min-forestcover-percent and forest-percent > 0 ]
  ask belowoptfc [ set do 3 ]
  let hr 0
  let do2patches patches with [ do = 2 ]
  ask do2patches [
    let aa count neighbors with [ do = 3 ]
    ifelse aa >= 1
    [ set hr forest-percent
      set hrng hrng + 1
      set dh hrng
      set do 1
      ask one-of neighbors with [ do = 3 ] [
        set hr hr + forest-percent
        set dh hrng
        set do 1
        ]
      ]
    [ set do 5 ]
    ]
  let do3patches patches with [ do = 3 ]
  ask do3patches [
  if any? neighbors with [ do = 3 ] [
    set hr forest-percent
    let nhr sum [ forest-percent ] of neighbors with [ do = 3 ]
    if (hr + nhr >= min-forestcover-percent) [
      set hrng hrng + 1
      set dh hrng
      set do 1
      ask neighbors with [ do = 3 ] [
        set dh hrng
        set do 1 ]
        ]
      ]
    ]
  ask patches with [ do > 1 ] [ set do -1 ]
  let tnhr 0
  set tnhr (max [ dh ] of patches)
  ;ask patches with [dh > 0] [ set pcolor scale-color blue dh 1 tnhr]          ;activate the following code to visualize deer homerange distribution in the landscape
  ask patches with [ dh > 0 ] [ set pcolor scale-color green forest-percent 1 0 ]
;patch procedure- initiate deer population
  let hrn-patches patches with [ dh > 0 ]
  let mp (precision (100 / (100 + sexratio))4) * 100
  let madp (precision (mp * adultprop) 2)
  let myp (precision (mp * yearlingprop) 2)
  let fadp (precision ((1 - mp) * adultprop) 2)
  let fyp (precision ((1 - mp) * yearlingprop) 2)
  let cds 0
  let xx 0
  ask hrn-patches [
    set cds cds + 1
    ifelse cds > tnhr
    [ stop ]
    [ set xx 0
      while [ xx < post_harvest_density ]             ;post-hunting (Jan) deer density per sq.mile
      [ let current-patch one-of patches with [ dh = cds ]
        ask current-patch [
          sprout-deers 1 [
            set shape "deer"
            set color orange
            set size 1.5
            set groid -1
            set aim 0
            ifelse random-float 100 < mp              ;male:female ratio
            [ set sex 1                                 ;1=male
              set gr -2
              ifelse random-float mp < madp
              [ set aim (20 + ((1 + random 14) * 12))
                set mgroid -1 ]
              [ ifelse (random-float (mp - madp) < myp)
                [ set aim 20
                  set mgroid -1 ]
                [ set aim 8
                  set mgroid -2 ]
                ]
              ]
            [ set sex 2                                   ;2=female
              ifelse random-float (1 - mp) < fadp
              [ set aim (20 + ((1 + random 14) * 12))
                ifelse random-float 1 < 0.5
                [ set gl 1 ]
                [ set gl 0 ]
                ]
              [ ifelse random-float ((1 - mp) - fadp) < fyp
                [ set aim 20 ]
                [ set aim 8 ]
                ]
              ]
            ht
            ]
          set xx xx + 1
          ]
        ]
      ]
    ]
  ;turtle procedure: doe social group formation
  let pot-leaders deers with [ gl = 1 and sex = 2 ]
  ask pot-leaders [
    set groid who
    set tgroid groid
    let pgr abs (round random-normal 2 1)
    let ndag count deers in-radius-nowrap 1.5 with [ sex = 2 and gl = 0 and gr = 0 ]
    if ndag > 0 [
      ifelse ndag >= pgr
      [ ask n-of pgr deers in-radius-nowrap 1.5 with [ sex = 2 and gl = 0 and gr = 0 ] [
        set gr -1
        set groid tgroid ]
        set gr pgr ]
      [ ask n-of ndag deers in-radius-nowrap 1.5 with [ sex = 2 and gl = 0 and gr = 0 ] [
        set gr -1
        set groid tgroid ]
        set gr ndag ]
      ]
    if gr < 2 [
      set gl 0
      set gr 0
      set groid -1
      ask deers with [ groid = tgroid ] [
        set gr 0
        set groid -1
        ]
      ]
    ]
  ask deers with [ sex = 2 and gl = 0 and groid = -1 and gr = 0 ] [ set gr -2 ]
 ;patch procedure: calculate mean forest percent for a patch and its immediate neighbors (required in the male yearling dispersal submodel)
  let nb-patches patches with [ border = 0 ]
  ask nb-patches [
    let tff forest-percent
    let nff sum [ forest-percent ] of neighbors with [ forest-percent > 0 ]
    set dfp (tff + nff) / 9 ]
  let b-patches patches with [ border >= 1 ]
  let vv 0
  ask b-patches [
    ask one-of neighbors with [ border = 0 ] [ set vv dfp ]
    set dfp vv
    ]
end
;-----------------------------------------------------------------------------------;
to-report min-forestcover-percent
  report 0.25
end
to-report max-forestcover-percent
  report 0.75
end
to-report doe-group-size-regulator
  report 6
end
to-report juvenile-pregnancy-rate
  report 20
end
to-report adult-pregnancy-rate
  report 80
end
to-report yearling-male-dispersal-rate                           ;dispersal rates for yearling males range between 46% to 80% (Long et al., 2005)
  report 0.46
end
to-report yearling-female-dispersal-rate
  report 0.22
end
to-report mean-female-dispersal-distance
 report 11
end
to-report stddev-dispersal-distance
 report 4
end
to-report mean-bachelor-group-size
 report 2 + random 4
end
to-report d
  report (remainder (ticks) 12) + 1
end
to-report year
  report (floor (ticks / 12) + 1)
end
to-report mad
  let ini-madults deers with [ sex = 1 and aim > 31 ]
  report count ini-madults
end
to-report fad
  let ini-fadults deers with [ sex = 2 and aim > 31 ]
  report count ini-fadults
end
to-report myd
  let ini-myearlings deers with [ sex = 1 and aim = 20 ]
  report count ini-myearlings
end
to-report fyd
  let ini-fyearlings deers with [ sex = 2 and aim = 20 ]
  report count ini-fyearlings
end
to-report mfd
  let ini-mfawns deers with [ sex = 1 and aim = 8 ]
  report count ini-mfawns
end
to-report ffd
  let ini-ffawns deers with [ sex = 2 and aim = 8 ]
  report count ini-ffawns
end
to-report mf
  let male-fawns deers with [ sex = 1 and aim < 12.5 ]
  report count male-fawns
end
to-report ff
  let female-fawns deers with [ sex = 2 and aim < 12.5 ]
  report count female-fawns
end
to-report my
  let male-yearlings deers with [ sex = 1 and aim > 12 and aim < 25 ]
  report count male-yearlings
end
to-report fy
  let female-yearlings deers with [ sex = 2 and aim > 12 and aim < 25 ]
  report count female-yearlings
end
to-report ma
  let male-adults deers with [ sex = 1 and aim > 24 ]
  report count male-adults
end
to-report fa
  let female-adults deers with [ sex = 2 and aim > 24 ]
  report count female-adults
end
;-----------------------------------------------------------------------------------;
to go
  set iteration (behaviorspace-run-number)
  if ticks = 312 [
    if output2 = "postharvest_population" or output2 = "both" [
      export-world (word "../results/PostHarvestPopulation" region "_v2.2.0.csv")
      ]
    stop
    ]
  if max [ gr ] of deers with [ sex = 2 and gl = 1 ] > 11 [
      print (word "min doe group size " min [ gr ] of deers with [ sex = 2 and gl = 1 ] )      ;to keep a check on group size during model run
      print (word "max doe group size " max [ gr ] of deers with [ sex = 2 and gl = 1 ] )
    ]
  ask deers [               ;turtle procedure: deer grow (age increases by 1 month) and deer die due to non-hunting causes of mortality
    st
    if random 100 < 90 [ ht ]
    individual-growth
    deer-die
  ]
  ask deers with [ sex = 2 and gl = 1 ] [
    set tgroid groid
    if gr < 2 [
      set gl 0
      set groid -1
      ask deers with [ groid = tgroid ] [
        set gr -2
        set groid -1
        set n_leaders_lost n_leaders_lost + 1
      ]
    ]
  ]
  if d > 1 and d < 10 [                                                ;turtle procedure: bachelor group leaders assess their group membership
    let male-leaders deers with [ sex = 1 and ml = 1 ]
    ask male-leaders [
      set tmgroid mgroid
      set gr count deers with [ sex = 1 and mgroid = tmgroid ]
      if gr <= 1 [ set ml 0 ]
    ]
    if min [ gr ] of deers with [ sex = 1 and ml = 1 ] < 2 or max [ gr ] of deers with [ sex = 1 and ml = 1 ] > 8 [
      print (word "min bachelor group size " min [ gr ] of deers with [ sex = 1 and ml = 1 ] )  ;to keep a check on group size during model run
      print (word "max bachelor group size " max [ gr ] of deers with [ sex = 1 and ml = 1 ] )
    ]
    ]
  if d = 10 [              ;turtle procedure: bachelor groups break down before the rutting season
    let male-leaders deers with [ sex = 1 and ml = 1 ]
    ask male-leaders [ set gr 0 ]
    ]
  if d = 1 [                                                          ;Post-harvest census
    let posthpop (mf + ff + my + fy + ma + fa)
    set vals1 (list (year) (mf) (my) (ma) (ff) (fy) (fa) (posthpop))
    let tmbg (my + ma) / 4
    let male-leaders deers with [ sex = 1 and ml = 1 ]                  ;turtle procedure: formation of bachelor groups
    let count-ml count male-leaders
    let pot-leaders deers with [ sex = 1 and aim > 32 and ml = 0 ]
    if tmbg > count-ml [
      ask n-of (tmbg - count-ml) pot-leaders [
        set ml 1
        set mgroid who ]
      ]
    set male-leaders deers with [ sex = 1 and ml = 1 ]
    ask male-leaders [ bachelor-group-formation ]
    ask male-leaders [
      set tmgroid mgroid
      set gr count deers with [ sex = 1 and mgroid = tmgroid ]
      ]
    let oldmales deers with [ sex = 1 and aim >= 229 ]
    let oldfemales deers with [ sex = 2 and aim >= 229 ]
    set oldm precision (count oldmales / ma) 3
    set oldf precision (count oldfemales / fa) 3
    ]
  if d = 5 [
    let male-yearlings deers with [ aim = 13 and sex = 1 ]
    ask male-yearlings [                                              ;turtle procedure: dispersal of male yearlings
      set tgroid groid
      set counter1 0
      if gr = -1 [ review-group-dynamics ]
      set gr -2
      set groid -1
      if random-float 1 < yearling-male-dispersal-rate [
        set mgroid -1
        deer-mdisperse
        ]
      ]
    let female-yearlings deers with [ aim = 13 and sex = 2 ]
    ask female-yearlings [                                            ;turtle procedure: dispersal of female yearlings
      if random-float 1 < yearling-female-dispersal-rate [
        set tgroid groid
        set counter1 0
        if gr = -1 [ review-group-dynamics ]
        set gr -2
        set groid -1
        deer-fdisperse
        ]
      ]
    let breeding-females deers with [ sex = 2 and aim > 12 ]
    ask breeding-females [
      let grfis 0
      if aim = 13 [
        if random 100 <= juvenile-pregnancy-rate + 1 [
          set grfis 0
          set tgroid groid
          ifelse tgroid >= 0
          [ ask deer tgroid [           ;turtle procedure: doe social group size adjustment in response to fawning
            ifelse gr > doe-group-size-regulator - 2
            [ set gr gr - 1
              set tgroid -1
              set tgr -2
              set grfis 1 ]
            [ set gr gr + 1
              set tgroid groid
              set tgr -1 ]
            ]
            if grfis = 1 [
              set groid -1
              set gr -2
              ]
            ]
          [ set tgroid -1
            set tgr -2 ]
          deer-reproduce                                      ;turtle procedure: female fawn breeding
          ]
        ]
      if aim > 24 [
        if random 100 < adult-pregnancy-rate + 1 [
          set grfis 0
          set tgroid groid
          ifelse gl > 0
          [ set gr gr + 2
            if gr > doe-group-size-regulator [               ;turtle procedure: new group formation as a response to doe breeding
              let xgr gr - doe-group-size-regulator
              let group-members deers with [ groid = tgroid and sex = 2 and aim > 13 and gl = 0 ]
              let ngr 0
              ifelse count group-members >= xgr
              [ set ngr xgr]
              [ set ngr count group-members ]
              ask n-of ngr group-members [ new-group-formation ]
              ]
            set tgr -1
            ]
          [ ifelse groid < 0
            [ ifelse aim > 36 and n_leaders_lost > 0  ;to prevent declining numbers of females with (gl = 1 )
              [ set gl 1
                set gr 2
                set groid who
                set tgroid who
                set tgr -1
                set n_leaders_lost n_leaders_lost - 1
                ]
              [ set tgroid -1
                set tgr -2 ]
              ]
            [ ask deer tgroid [
                ifelse gr > 4
                [ set gr gr - 1
                  set tgroid -1
                  set tgr -2
                  set grfis 1
                ]
                [ set gr gr + 2
                  set tgroid groid
                  set tgr -1
                ]
              ]
              if grfis = 1 [
                set groid -1
                set gr -2
                ]
              ]
            ]
          deer-reproduce                                         ;turtle procedure: adult doe breeding
          ]
        ]
                                                                 ;turtle procedure: adjustment of doe social group size after the breeding season
      if gl = 1 and gr < 4 [
        set tgroid groid
        let solitary-adult-females-here deers in-radius-nowrap 1.5 with [ sex = 2 and gr = -2 and aim >= 13 ]
        let sd count solitary-adult-females-here
        let sd1 0
        set tgr 0
        if sd > 0 [
          ifelse sd > 2
          [ set sd1 2 ]
          [ set sd1 1 ]
          ask n-of sd1 solitary-adult-females-here [
            let tmomid who
            set groid tgroid
            set gr -1
            set tgr tgr + 1
            if any? deers with [ momid = tmomid and aim = 1 ] [
              let my-fawns deers with [ momid = tmomid and aim = 1 ]
              ask my-fawns [
                set groid tgroid
                set gr -1
                set tgr tgr + 1
                ]
              ]
            ]
          ]
        set gr gr + tgr
        if gr = 0 [
          set gl 0
          set groid -1
          set gr -2
          ]
        ]
      ]
    ]
  if d = 11 [
    let solitary-male-yearlings deers with [ sex = 1 and aim = 19 and gr = -2 and mgroid = -2 ]
    ask solitary-male-yearlings [                            ;turtle procedure: male yearling dispersal before the rutting season
      if random-float 1 < yearling-male-dispersal-rate [
        set mgroid -1
        deer-mdisperse ]
      set mgroid -1
      ]
    if year >= 2 [                                           ;for lambda calculation
      set nt ntplus1
      ]
    let phn (mf + my + ma + ff + fy + fa)                    ;Pre-harvest census
    set ntplus1 phn                                          ;for lambda calculation
    ifelse year >= 2
    [ set lambda precision (ntplus1 / nt) 3
    ]
    [ set lambda "NA" ]
    set vals2 (list (mf) (my) (ma) (ff) (fy) (fa) (phn) (lambda))
    if ticks = 310 [
      if (output2 = "preharvest_population" or output2 = "both") [
        export-world (word "../results/PreHarvestPopulation" region "_v2.2.0.csv")
        ]
      ]
      ;sensitivity analysis
    file-open (word "../results/sa" region ".csv")
    file-type year
    file-type ","
    file-type phn
    file-type ","
    file-type precision ((ma + fa) / phn) 2
    file-type ","
    file-type precision ((my + fy) / phn) 2
    file-type ","
    file-type precision ((mf + ff) / phn) 2
    file-type ","
    file-type precision ((ff + fy + fa) / (mf + my + ma)) 2
    file-type ","
    file-type lambda
    file-type ","
    file-type iteration
    file-print ""
    file-close
    ]
  if d = 12 [
    ask deers [                                             ;turtle procedure: harvest mortality
      if aim < 10 [
        ifelse sex = 1
        [ if random-float 1 < mf12hm [
          set tgroid groid
          hunting-mortality-mf12
          ]
          ]
        [ if random-float 1 < ff12hm [
          set tgroid groid
          set twho who
          hunting-mortality-ff12
          ]
          ]
        ]
      if aim = 20 [
        ifelse sex = 1
        [ if random-float 1 < myhm [
          set tgroid groid
          hunting-mortality-my
          ]
          ]
        [ if random-float 1 < fyhm [
          set tgroid groid
          set twho who
          hunting-mortality-fy
          ]
          ]
        ]
      if aim > 30 [
        ifelse sex = 1
        [ if random-float 1 < mahm [
          set tgroid groid
          hunting-mortality-ma
          ]
          ]
        [ if random-float 1 < fahm [
          set tgroid groid
          set twho who
          hunting-mortality-fa
          ]
          ]
        ]
      ]
    let tot_harvest (tmfh + tmyh + tamh + tffh + tfyh + tafh)
    set vals3 (list (tmfh) (tmyh) (tamh) (tffh) (tfyh) (tafh) (tot_harvest) (iteration))
    set vals (sentence vals1 vals2 vals3)
    file-open (word "../results/deerpopdy" region ".csv")
    file-type first vals
    foreach but-first vals [ [ ?1 ] ->
      file-type "," file-type ?1
      ]
    file-print""
    file-close
    set ndd 0
    if ticks < 310 [
      set tamh 0
      set tafh 0
      set tmfh 0
      set tffh 0
      set tmyh 0
      set tfyh 0
      ]
    ]
  set-current-plot "deer population"
  plotxy ticks count deers
  tick
end
;turtle procedure: age (age in months) increases by 1 time step.
to individual-growth
  set aim aim + 1
end
;turtle procedure: female deer reproduce
to deer-reproduce
  let mom who
  ifelse aim < 13.5
  [ hatch-deers 1 [
    set aim 1
    set shape "deer"
    set color orange
    set size 1.5
    set gl 0
    set momid mom
    set groid tgroid
    set gr tgr
    ifelse random 100 < 51
    [ set sex 1
      set mgroid -2 ]
    [ set sex 2 ]
    ]
  ]
  [ hatch-deers 2 [
    set aim 1
    set shape "deer"
    set color orange
    set size 1.5
    set gl 0
    set momid mom
    set groid tgroid
    set gr tgr
    ifelse random 100 < 51
    [ set sex 1
      set mgroid -2 ]
    [ set sex 2 ]
    ]
  ]
end
to deer-mdisperse                                         ;turtle procedure: male yearling dispersal
  ask patch-here [
    if dfp < 0.72 [
      let mdd (35.07 - (48.14 * (dfp)))                   ;relationship between mean dispersal distance and proportion of forest cover to estimate the mean dispersal distance Long et al., 2005; Difenbach et al., 2008
      let sddd sqrt(e ^ (3.51 + (0.077 * mdd)))
      let lv ln (1 + (sddd ^ 2) / (mdd ^ 2))
      let lm ln mdd - (lv / 2)
      let ls sqrt lv
      set dd (exp (random-normal lm ls) * 0.6214) ]       ;this is in miles, 1 patch fd is 1 mile, so fd dd
    ]
  let counter-md 0
  rt random 360
  while [ counter-md < round dd ]
  [ ask patch-here [
      set dxp (pxcor)
      set dyp (pycor) ]
    fd 1
    set counter-md counter-md + 1
    ask patch-here [
      set dxn (pxcor)
      set dyn (pycor) ]
    if (abs(dxn - dxp) > 1 or abs(dyn - dyp) > 1)
    [ set ndd ndd + 1
      set momid 0 ]
  ]
  finalize-home-patch
end
to deer-fdisperse                                         ;turtle procedure: female yearling dispersal
  let counter-fd 0
  rt random 360
  set dd round (random-normal mean-female-dispersal-distance stddev-dispersal-distance)
  while [ counter-fd < dd ]
  [ ask patch-here [
      set dxp pxcor
      set dyp pycor ]
    fd 1
    set counter-fd counter-fd + 1
    ask patch-here [
      set dxn pxcor
      set dyn pycor ]
    if (abs(dxn - dxp) > 1 or abs(dyn - dyp) > 1)
    [ set ndd ndd + 1
      set momid 0 ]
    ]
  finalize-home-patch
end
to finalize-home-patch                                   ;turtle procedure: if a deer ends up on a non-deer habitat patch after dispersing, it goes to the nearest deer habitat patch
  let fhp 0
  ask patch-here [
    if do != 1 [
      set fhp 1 ]
    ]
  if fhp > 0 [
    set dxp pxcor
    set dyp pycor
    move-to min-one-of patches with [ do = 1 ] [ distance myself ]
    set color blue
    ask patch-here [
      set dxn pxcor
      set dyn pycor ]
    if (abs(dxn - dxp) > 1 or abs(dyn - dyp) > 1)        ;dispersing deer goes out of the model landscape
    [ set ndd ndd + 1
      set momid 0 ]
    ]
end
to new-group-formation                                   ;turtle procedure: new group formation after breeding season
  set tgroid groid
  set groid -1
  set ttgroid who
  set tgr 0
  if any? deers with [ momid = ttgroid and aim = 1 ] [
    let my-fawns deers with [ momid = ttgroid and aim = 1 ]
    ask my-fawns [
      set gr -2
      set groid -1
      set tgr tgr + 1 ]
    ]
  set gr -2
  ask deer tgroid [ set gr gr - (tgr + 1) ]
end
to deer-die                                            ;turtle procedure: deer non-harvest mortality
  set tgroid groid
  if sex = 2 [ set twho who ]
  ;fawns upto 6 months
  ifelse aim < 6.5
  [ ifelse sex = 1
    [ if precision (random-float 1) 3 < mf6nhm [
      set counter1 0
      if gr = -1 [ review-group-dynamics ]
      die
      ]
      ]
  [ if precision (random-float 1) 3 < ff6nhm [
      set counter1 0
      if gr = -1 [ review-group-dynamics ]
      die
      ]
      ]
    ]
  ;fawns 7 to 12 months
  [ ifelse aim < 12.5
    [ ifelse sex = 1
      [ if random-float 1 < mf12nhm [
        set counter1 0
        if gr = -1 [ review-group-dynamics ]
        die
        ]
        ]
      [ if random-float 1 < ff12nhm [
        set counter1 0
        if gr = -1 [ review-group-dynamics ]
        die
        ]
        ]
      ]
    ;yearlings 13 to 24 months
    [ ifelse aim < 24.5
      [ ifelse sex = 1
        [ if random-float 1 < mynhm [
          set counter1 0
          if gr = -1 [ review-group-dynamics ]
          die
          ]
          ]
        [ if random-float 1 < fynhm [
          if any? deers with [ momid = twho and aim < 2.5 ] [      ;turtle procedure: fawns less than 2 months old die if their mother dies
            let my-fawns deers with [ momid = twho and aim < 2.5 ]
            ask my-fawns [
              set counter1 counter1 + 1
              die
              ]
            ]
          ifelse gl > 0
          [ let pot-groupleaders deers with [ groid = tgroid and gl = 0 and sex = 2 and aim >= 18 ] ;turtle procedure: if a doe social group leader dies, leadership is transferred or the group breaks down
            ifelse count pot-groupleaders > 0
            [ ask one-of pot-groupleaders [
              set ttgroid who
              set gl 1
              set groid ttgroid
              let transfer-group deers with [ groid = tgroid and gl = 0 ]
              if count transfer-group > 0 [
                ask transfer-group [
                  set groid ttgroid
                  ]
                ]
              let new-group deers with [ groid = ttgroid and gl = 0 ]
              set gr count new-group
              ]
              ]
            [ let other-groupleaders-here deers-here with [ gl = 1 and groid != tgroid and gr < 3 ]
              ifelse count other-groupleaders-here > 0
              [ ask one-of other-groupleaders-here [
                set ttgroid groid
                let transfer-deers deers with [ groid = tgroid and gl = 0 ]
                set gr gr + count transfer-deers
                ]
                let new-group deers with [ groid = tgroid and gl = 0 ]
                ask new-group [ set groid ttgroid ]
                ]
              [ let group-members deers with [ groid = tgroid and gl = 0 ]
                ask group-members [
                  set gr -2
                  set groid -1
                  ]
                ]
              set n_leaders_lost n_leaders_lost + 1
              ]
            ]
          [ if gr = -1 [ review-group-dynamics ]
            ]
          set counter1 0
          die
          ]
          ]
        ]
      ;male 25 to 240 and more than 240
      [ ifelse sex = 1
        [ ifelse aim < 240
          [ if random-float 1 < precision (manhm - oldm) 3 [
            if ml = 1 [
              set tmgroid mgroid
              let my-group deers with [ mgroid = tmgroid ]
              ask my-group [
                set mgroid -1
                ]
              ]
            die
            ]
            ]
          [ if random-float 1 < 0.8 [
            if ml = 1 [
              set tmgroid mgroid
              let my-group deers with [ mgroid = tmgroid ]
              ask my-group [
                set mgroid -1
                ]
              ]
            die
            ]
            ]
          ]
        ;female 25 to 240
        [ ifelse aim < 240
          [ if random-float 1 < precision (fanhm - oldf) 3 [
            let my-fawns deers with [ momid = twho and aim < 2.5 ]
            if count my-fawns > 0 [                               ;turtle procedure: fawns less than 2 months old die if their mother dies
              ask my-fawns [
                set counter1 counter1 + 1
                die
                ]
              ]
            ifelse gl > 0
            [ let pot-groupleaders deers with [ groid = tgroid and aim > 18 and sex = 2 and gl = 0 ]  ;turtle procedure: if a doe social group leader dies, leadership is transferred or the group breaks down
              let zz count pot-groupleaders
              ifelse zz > 0
              [ ifelse any? pot-groupleaders with [ aim > 29 ]
                [ ask one-of pot-groupleaders with [ aim > 29 ] [
                  set gl 1
                  set ttgroid who
                  set groid ttgroid
                  let my-group deers with [ groid = tgroid and gl = 0 ]
                  ask my-group [ set groid ttgroid ]
                  let new-group deers with [ groid = ttgroid and gl = 0 ]
                  set gr count new-group
                  ]
                  ]
                [ ask one-of pot-groupleaders [
                  set gl 1
                  set ttgroid who
                  set groid ttgroid
                  let my-group deers with [ groid = tgroid and gl = 0 ]
                  ask my-group [ set groid ttgroid ]
                  let new-group deers with [ groid = ttgroid and gl = 0 ]
                  set gr count new-group
                  ]
                  ]
                ]
              [ let other-groupleaders-here deers-here with [ gl = 1 and groid != tgroid and gr < 3 ]
                ifelse count other-groupleaders-here > 0
                [ ask one-of other-groupleaders-here [
                  set ttgroid groid
                  let add-deers deers with [ groid = tgroid and gl = 0 ]
                  set gr gr + count add-deers
                  ]
                  let new-group deers with [ groid = tgroid and gl = 0 ]
                  ask new-group [ set groid ttgroid ]
                  ]
                [ let group-members deers with [ groid = tgroid and gl = 0 ]
                  ask group-members [
                    set gr -2
                    set groid -1
                    ]
                  ]
                set n_leaders_lost n_leaders_lost + 1
                ]
              ]
            [ if gr = -1 [ review-group-dynamics ] ]
            set counter1 0
            die
            ]
            ]
          ;female 240 and more
          [ if random-float 1 < 0.8 [
            if any? deers with [ momid = twho and aim < 2.5 ] [               ;turtle procedure: fawns less than 2 months old die if their mother dies
              let my-fawns deers with [ momid = twho and aim < 2.5 ]
              ask my-fawns [
                set counter1 counter1 + 1
                die
                ]
              ]
            ifelse gl > 0
            [ let pot-groupleaders deers with [ groid = tgroid and aim > 18 and sex = 2 and gl = 0 ]    ;turtle procedure: if a doe social group leader dies, leadership is transferred or the group breaks down
              let zz count pot-groupleaders
              ifelse zz > 0
              [ ifelse any? pot-groupleaders with [ aim > 29 ]
                [ ask one-of pot-groupleaders with [ aim > 29 ] [
                  set gl 1
                  set ttgroid who
                  set groid ttgroid
                  let my-group deers with [ groid = tgroid and gl = 0 ]
                  ask my-group[ set groid ttgroid ]
                  let new-group deers with [ groid = ttgroid and gl = 0 ]
                  set gr count new-group
                  ]
                  ]
                [ ask one-of pot-groupleaders [
                  set gl 1
                  set ttgroid who
                  set groid ttgroid
                  let my-group deers with [ groid = tgroid and gl = 0 ]
                  ask my-group [ set groid ttgroid ]
                  let new-group deers with [ groid = ttgroid and gl = 0 ]
                  set gr count new-group
                  ]
                  ]
                ]
              [ let other-groupleaders-here deers-here with [ gl = 1 and groid != tgroid and gr < 3 ]
                ifelse count other-groupleaders-here > 0
                [ ask one-of other-groupleaders-here [
                  set ttgroid groid
                  let my-group deers with [ groid = tgroid and gl = 0 ]
                  set gr gr + count my-group
                  ]
                  let new-group deers with [ groid = tgroid and gl = 0 ]
                  ask new-group [ set groid ttgroid ]
                  ]
                [ let my-group deers with [ groid = tgroid and gl = 0 ]
                  ask my-group [
                    set gr -2
                    set groid -1 ]
                  ]
                set n_leaders_lost n_leaders_lost + 1
                ]
              ]
            [ if gr = -1 [ review-group-dynamics ] ]
            set counter1 0
            die
            ]
            ]
          ]
        ]
      ]
    ]
end
;turtle procedure: male fawn harvest mortality
to hunting-mortality-mf12
  if gr = -1 [
    set counter1 0
    review-group-dynamics
    ]
  set tmfh tmfh + 1
  die
end
;turtle procedure: female fawn harvest mortality
to hunting-mortality-ff12
  if gr = -1 [
    set counter1 0
    review-group-dynamics
    ]
  set tffh tffh + 1
  die
end
;turtle procedure: male yearling harvest mortality
to hunting-mortality-my
  set tmyh tmyh + 1
  die
end
;turtle procedure: female yearling harvest mortality
to hunting-mortality-fy
  ifelse gl > 0
  [ let pot-groupleaders deers with [ groid = tgroid and gl = 0 and sex = 2 and aim >= 18 ]       ;turtle procedure: if a doe social group leader dies, leadership is transferred or the group breaks down
    ifelse count pot-groupleaders > 0
    [ ask one-of pot-groupleaders [
      set ttgroid who
      set gl 1
      set groid ttgroid
      let transfer-group deers with [ groid = tgroid and gl = 0 ]
      let xgr count transfer-group
      ask transfer-group [ set groid ttgroid ]
      set gr xgr
      ]
      ]
    [ let other-groupleaders-here deers with [ gl = 1 and groid != tgroid and gr < 3 ]
      ifelse count other-groupleaders-here > 0
      [ ask one-of other-groupleaders-here [
        set ttgroid groid
        let transfer-group deers with [ groid = tgroid and gl = 0 ]
        set gr gr + count transfer-group
        ]
        let new-group deers with [ groid = tgroid and gl = 0 ]
        ask new-group [ set groid ttgroid ]
        ]
      [ let group-members deers with [ groid = tgroid and gl = 0 ]
        ask group-members [
          set gr -2
          set groid -1
          ]
        ]
      set n_leaders_lost n_leaders_lost + 1
      ]
    ]
  [ if gr = -1 [
    set counter1 0
    review-group-dynamics
    ]
    ]
  set tfyh tfyh + 1
  die
end
;turtle procedure: adult male harvest mortality
to hunting-mortality-ma
  if ml = 1 [
    set tmgroid mgroid
    let group-members deers with [ mgroid = tmgroid and ml = 0 ]
    ask group-members [ set mgroid -1 ]
    ]
  set tamh tamh + 1
  die
end
;turtle procedure: adult female harvest mortality
to hunting-mortality-fa
  set tafh tafh + 1
  if gl = 1 [
    let pot-groupleaders deers with [ groid = tgroid and sex = 2 and aim > 18 and gl = 0 ]  ;turtle procedure: if a doe social group leader dies, leadership is transferred or the group breaks down
    let zz count pot-groupleaders
    ifelse zz > 0
    [ ifelse any? pot-groupleaders with [ aim > 29 ]
      [ ask one-of pot-groupleaders with [ aim > 29 ] [
        set gl 1
        set ttgroid who
        set groid ttgroid
        let transfer-group deers with [ groid = tgroid and gl = 0 ]
        let xtr count transfer-group
        ask transfer-group [ set groid ttgroid ]
        set gr xtr
       ]
      ]
      [ ask one-of pot-groupleaders [
        set gl 1
        set ttgroid who
        set groid ttgroid
        let transfer-group deers with [ groid = tgroid and gl = 0 ]
        let xtr count transfer-group
        ask transfer-group [ set groid ttgroid ]
        set gr xtr
        ]
        ]
      ]
    [ let other-groupleaders-here deers-here with [ gl = 1 and groid != tgroid and gr < 3 ]
      ifelse count other-groupleaders-here > 0
      [ ask one-of other-groupleaders-here [
       set ttgroid groid
       let transfer-group deers with [ groid = tgroid and gl = 0 ]
       let xtr count transfer-group
       set gr gr + xtr
       ask transfer-group [ set groid ttgroid ]
        ]
        ]
      [ let transfer-group deers with [ groid = tgroid and gl = 0 ]
        ask transfer-group [
          set gr -2
          set groid -1
        ]
        ]
      set n_leaders_lost n_leaders_lost + 1
      ]
    ]
  if gr = -1 [
    set counter1 0
    review-group-dynamics
    ]
  die
end
to bachelor-group-formation             ;turtle procedure: bachelor group formation
  set tmgroid mgroid
  let group-members deers in-radius-nowrap 1.5 with [ sex = 1 and mgroid = tmgroid and ml = 0 ]
  let cgs count group-members
  let cgs1 0
  let pot-gr-size mean-bachelor-group-size
  if cgs < pot-gr-size [
    let pot-groupmembers-here deers in-radius-nowrap 1.5 with [ sex = 1 and mgroid = -1 and ml = 0 ]
    set cgs1 count pot-groupmembers-here
    ifelse cgs1 >= (pot-gr-size - cgs)
    [ ask n-of (pot-gr-size - cgs) pot-groupmembers-here [
      set mgroid tmgroid
      ]
      ]
    [ if cgs1 > 0 [
      ask n-of cgs1 pot-groupmembers-here [
        set mgroid tmgroid
        ]
      ]
      ]
    ]
end
to review-group-dynamics              ;turtle procedure: doe social group leader loses leadership status if no group members left
  ask deer tgroid [
    set gr gr - (counter1 + 1)
    if gr <= 0 [
      set gl 0
      set groid -1
      set gr -2
      set n_leaders_lost n_leaders_lost + 1
      ]
    ]
end
@#$#@#$#@
GRAPHICS-WINDOW
798
65
1298
510
-1
-1
4.0
1
10
1
1
1
0
1
1
1
0
122
0
108
1
1
1
ticks
30.0

BUTTON
441
166
505
199
Setup
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

MONITOR
720
65
787
110
total deer
count deers
17
1
11

SLIDER
158
458
277
491
mf6nhm
mf6nhm
0
.1
0.055
.001
1
NIL
HORIZONTAL

SLIDER
159
497
278
530
ff6nhm
ff6nhm
0
.1
0.055
.001
1
NIL
HORIZONTAL

SLIDER
158
535
278
568
mf12nhm
mf12nhm
0
1
0.05
.01
1
NIL
HORIZONTAL

SLIDER
157
575
280
608
ff12nhm
ff12nhm
0
1
0.05
.01
1
NIL
HORIZONTAL

SLIDER
303
617
424
650
myhm
myhm
0
1
0.25
.01
1
NIL
HORIZONTAL

SLIDER
303
656
424
689
fyhm
fyhm
0
1
0.15
.01
1
NIL
HORIZONTAL

SLIDER
157
614
278
647
mynhm
mynhm
0
1
0.01
.01
1
NIL
HORIZONTAL

SLIDER
156
653
278
686
fynhm
fynhm
0
1
0.0
0.01
1
NIL
HORIZONTAL

SLIDER
302
693
425
726
mahm
mahm
0
1
0.4
.01
1
NIL
HORIZONTAL

SLIDER
303
732
426
765
fahm
fahm
0
1
0.2
.001
1
NIL
HORIZONTAL

SLIDER
156
691
279
724
manhm
manhm
0
1
0.01
.01
1
NIL
HORIZONTAL

SLIDER
155
731
279
764
fanhm
fanhm
0
1
0.02
.01
1
NIL
HORIZONTAL

BUTTON
442
235
505
268
Go
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

PLOT
626
113
786
233
deer population
months
deer
0.0
10.0
0.0
15000.0
true
false
"" ""
PENS
"pen-0" 1.0 0 -2674135 true "" ""

SLIDER
301
458
422
491
mf6hm
mf6hm
0
0
0.0
0
1
NIL
HORIZONTAL

SLIDER
302
500
423
533
ff6hm
ff6hm
0
0
0.0
0
1
NIL
HORIZONTAL

SLIDER
302
539
424
572
mf12hm
mf12hm
0
1
0.05
.01
1
NIL
HORIZONTAL

SLIDER
303
577
423
610
ff12hm
ff12hm
0
1
0.02
.01
1
NIL
HORIZONTAL

PLOT
629
364
789
484
Doe group size
NIL
NIL
2.0
15.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 1 -16777216 true "" "histogram [gr + 1] of deers with [gl = 1]"

PLOT
626
237
786
357
Bachelor group size
NIL
NIL
0.0
15.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 1 -16777216 true "" "histogram [ gr ] of deers with [ ml > 0 ]"

TEXTBOX
40
460
148
492
Young male fawns
13
105.0
1

TEXTBOX
160
421
293
453
Non-hunting mortality\n(monthly rates)
13
105.0
1

TEXTBOX
306
422
415
454
Hunting mortality\n(annual rates)
13
105.0
1

TEXTBOX
32
499
151
531
Young female fawns
13
105.0
1

TEXTBOX
44
537
151
569
Older male fawns
13
105.0
1

TEXTBOX
31
578
146
610
Older female fawns
13
105.0
1

TEXTBOX
60
616
147
648
Male yearlings
13
105.0
1

TEXTBOX
45
654
148
686
Female yearlings
13
105.0
1

TEXTBOX
76
692
143
710
Male adults
13
105.0
1

TEXTBOX
61
733
145
765
Female adults
13
105.0
1

CHOOSER
180
85
391
130
region
region
"Elkhart" "LaGrange" "LaPorte" "StJoseph" "Kankakee"
4

SLIDER
179
238
391
271
post_harvest_density
post_harvest_density
0
50
14.0
1
1
per sq. mile
HORIZONTAL

SLIDER
179
284
392
317
sexratio
sexratio
0
200
160.0
10
1
females per 100 males
HORIZONTAL

SLIDER
178
327
392
360
adultprop
adultprop
0
1
0.45
.05
1
NIL
HORIZONTAL

SLIDER
177
372
392
405
yearlingprop
yearlingprop
0
1
0.25
.05
1
NIL
HORIZONTAL

SWITCH
181
164
390
197
recommended_parameter_values
recommended_parameter_values
1
1
-1000

CHOOSER
437
85
550
130
output2
output2
"preharvest_population" "postharvest_population" "both"
2

MONITOR
977
16
1034
61
Month
d
17
1
11

MONITOR
918
16
968
61
Year
year
17
1
11

MONITOR
800
15
910
60
NIL
Region
17
1
11

TEXTBOX
184
201
382
226
--------- or -----------
22
0.0
1

TEXTBOX
180
54
330
73
1. Select the region.
15
15.0
1

TEXTBOX
180
137
376
160
2. Set parameters.
15
15.0
1

TEXTBOX
436
57
586
76
3. Select output files.
15
15.0
1

TEXTBOX
443
140
508
160
4. Setup.
15
15.0
1

TEXTBOX
444
210
566
229
5. Run the model.
15
15.0
1

TEXTBOX
10
151
160
303
If 'On', parameter values extrapolated for Midwest deer populations will be used. Alternatively, sliders can be used to set parameter values.
15
15.0
1

TEXTBOX
5
10
610
33
INdiana Odocoileus virginianus POPulation (INOvPOP) ABM version 2.2.0
18
0.0
1

@#$#@#$#@
## WHAT IS IT?

**MOOvPOP** simulates population dynamics of white-tailed deer and generates pre-harvest deer population (abundance, sex-age composition and distribution in the landscape) for a selected sampling region.This model generated pre-harvest population can be used to initialize another model **MOOvPOPsurveillance** which simulates hunter harvest and CWD testing under different assumptions (random CWD distribution in the deer population/clustered CWD distribution in the deer population; random sampling / non-random sampling).

## HOW IT WORKS

Processes like social organization, group dynamics, dispersal, and hunting mortality occur at an individual level and influence interactions among individuals. Such interactions underpin host heterogeneity, and thereby influence disease transmission in a host population. We have incorporated these processes in **MOOvPOP** so that the model-generated population reflects heterogeneity observed in real-world host populations. 

## HOW TO USE IT

User selects a sampling region (Chooser: _Region_), after making sure that the GIS file for the desired sampling region is available in the same folder as that of the model. GIS files for some Missouri Counties are available [here] (https://github.com/anyadoc/MOOvPOP). Use the ‘Clone or download’ tab and remember to unzip the folder after downloading.  These files should be in the **same** folder as the NetLogo model file (**MOOvPOP**).

User then sets the values for different parameters using the sliders provided on the interface ( _PostHarvestDensity_ , _sexratio_ , _adultprop_ , _yearlingprop_ , and gender- and age-class wise _hunting_ and _non-hunting_ mortality rates). Alternatively, if the Switch _Recommended_parameter_values_ is on, the model will use pre-determined parameter values for the selected region. 

User should provide the file path and name in the input window _'PopulationDynamicsOutputFile'_. This output file (.csv) documents annual deer population dynamics (gender- and age-class wise pre-harvest,post-harvest deer abundance, as well as harvest distribution). The user further selects desired outputs (Chooser:  _Ouput2_ ): _preharvest_population_ , _postharvest_population_ or _both_ . Based on the selection, in the last year of model run, respective deer populations are stored as .csv files (export-world).


## THINGS TO NOTICE

The interface has several monitors and graphical outputs for assessing the model performance. Six monitors ( _male adult deer_ , _male yearling deer_ , _male fawn deer_ , _female adult deer_ , _female yearling deer_ and _female fawn deer_ ) document the initial gender- and age-specific composition of the deer population. The graph 'deer population' plots deer abundance at every time step, and two histograms 'bachelor group size' and 'doe group size' present the group size distribution in the deer population. These plots and monitors provide a mechanism to assess the functioning of different model processess like harvest mortality and grouping behavior of white-tailed deer.


## THINGS TO TRY

Change hunting and/or non-hunting mortality of different age-/gender-classes and notice the effect on overall deer abundance. 

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

[MOOvPOPsurveillance] (https://www.openabm.org/model/5576/version/4/view)

## CREDITS AND REFERENCES

Belsare, Aniruddha V., Gompper, Matthew E., Millspaugh, Joshua J. (2017, July 24). "MOOvPOP" (Version 5). CoMSES Computational Model Library. Retrieved from: https://www.openabm.org/model/5585/version/5 
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

deer
false
0
Polygon -7500403 true true 195 210 210 255 195 240 180 195 165 165 135 165 105 165 75 165 72 211 60 210 60 180 45 150 45 120 30 90 45 105 180 105 225 45 225 60 270 90 255 90 225 90 180 150
Polygon -7500403 true true 73 210 86 251 75 240 60 210
Polygon -7500403 true true 45 105 30 75 30 90 45 105 60 120 45 120
Line -7500403 true 210 60 165 15
Line -7500403 true 225 60 255 45
Line -7500403 true 195 45 210 15
Line -7500403 true 255 45 255 30
Line -7500403 true 255 45 270 30
Line -7500403 true 195 15 180 30

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
NetLogo 6.2.2
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="calibrationLaporte" repetitions="10" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="recommended_parameter_values">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="manhm">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mf6nhm">
      <value value="0.055"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sexratio">
      <value value="160"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ff12hm">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="yearlingprop">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ff6nhm">
      <value value="0.055"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="output2">
      <value value="&quot;both&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ff12nhm">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mf12nhm">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fahm">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fynhm">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="region">
      <value value="&quot;LaPorte&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="adultprop">
      <value value="0.45"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mahm">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ff6hm">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="myhm">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fanhm">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fyhm">
      <value value="0.15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mf12hm">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mynhm">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="post_harvest_density">
      <value value="14"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mf6hm">
      <value value="0"/>
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
