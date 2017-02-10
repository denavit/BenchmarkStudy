set startTime [clock seconds]
set testOutputFlag 1

if {$frameType == "Sidesway_Uninhibited"} {

} elseif {$frameType == "Sidesway_Inhibited"} {
	set Delta0 0
	set gamma  0
} else {
    puts "frameType: $frameType not recgonized"
}

# ############ BUILD MODEL ############
# Define nodes and constraints for column
set PI [expr 2*asin(1.0)]
for { set i 0 } { $i <= $numEles } { incr i 1 } {
    set x [expr double($i)/double($numEles)]
    set imperf [expr $Delta0*$x + $delta0*sin($PI*$x)]
    node [expr $i+1] $imperf [expr $x*$Lc]
}


fix 1 1 1 0

if {$frameType == "Sidesway_Uninhibited"} {

	# Define nodes and constraints for rotational springs
	if { $rotStiffBot != "Fixed" && $rotStiffBot != "Free" } {
	    node 200 0 0
	    fix 200 1 1 1 
	} else {
    	if { $rotStiffBot == "Fixed" } {
	        fix 1 0 0 1
    	}
	}

	if { $rotStiffTop != "Fixed" && $rotStiffTop != "Free" } {
    	set imperf $Delta0    
	    node 201 $imperf $Lc
    	fix 201 1 1 1 
	} else {
	    if { $rotStiffTop == "Fixed" } {
    	    fix $numNodes 0 0 1
	    }
	}

	# Define elements for rotational springs
	if { $rotStiffBot != "Fixed" && $rotStiffBot != "Free" } {
	    uniaxialMaterial Elastic 200 $rotStiffBot
	    element zeroLength 200 1 200 -mat 200 -dir 6
	} 
	if { $rotStiffTop != "Fixed" && $rotStiffTop != "Free" } {
	    uniaxialMaterial Elastic 201 $rotStiffTop
    	element zeroLength 201 $numNodes 201 -mat 201 -dir 6 
	} 


} elseif {$frameType == "Sidesway_Inhibited"} {
	fix $numNodes 1 0 0
	
} else {
    puts "frameType: $frameType not recgonized"
}

# Define Geometric Transformation
set geomTransfTag 1
geomTransf $geomTransfType $geomTransfTag

# Define elements for column 
set extraArgs [list]
if {$tryNumber >= 2} {
  # lappend extraArgs -geomLinear
}
for { set i 1 } { $i <= $numEles } { incr i 1 } {
  element mixedBeamColumn2d $i $i [expr $i+1] $numIP $columnSectionTag $geomTransfTag {*}$extraArgs
}

# Define nodes and constraints for leaning columns
if { $gamma != 0 } {
   	set imperf $Delta0
    node 100 0 0
    node 101 $imperf $Lc
    fix 100 1 1 0
    equalDOF $numNodes 101 1
    element elasticBeamColumn 100 100 101 $leaningA $leaningE $leaningI $geomTransfTag
}

# ############ DEFINE RECORDERS ############
recorder Node    -file $nodeDisplacementFilename -nodeRange 1 $numNodes -dof 1 2 3 disp
recorder Element -file $elementForceFilename -eleRange 1 $numEles localForce
recorder Element -file $elementSectionDeformationFilename -eleRange 1 $numEles sectionDeformation_Force
recorder Element -file $elementSectionStiffnessFilename -eleRange 1 $numEles sectionStiffness
set eigenFileId [open $eigenFilename w]


# ############ LOADING AND ANALYSIS OPTIONS ############
# Analysis Options
system UmfPack  
constraints Transformation
test NormUnbalance $baseForceTolerance 20 $testOutputFlag
algorithm Newton 
numberer Plain

# One step with nothing
integrator LoadControl 0
analysis Static
set ok [analyze 1] 
if {$ok == 0} {
    set lowestEigenValue [eigen -standard -symmBandLapack 1]
    set firstEigenValue $lowestEigenValue
    puts $eigenFileId "[getTime] $lowestEigenValue" 
}

# ############ DISPLACEMENT CONTROL / AXIAL ONLY LOADING ############
if {$analysisType == "LimitPoint_Proportional"} {
    # Apply gravity loads
    pattern Plain 1 Linear {
    	if {$frameType == "Sidesway_Uninhibited"} {
            load $numNodes $H -1.0 0.0  
            if { $gamma != 0 } { 
                load  101 0.0 [expr -1*$gamma] 0.0  
            }
		} elseif {$frameType == "Sidesway_Inhibited"} {         
            load 1         0.0  0.0 [expr  -1.0*$M]
        	load $numNodes 0.0 -1.0 [expr $beta*$M]
		} else {
		    puts "frameType: $frameType not recgonized"
		}
    }

    set ok 0
    set iStep 0


	if {$frameType == "Sidesway_Uninhibited"} {
		set controlledNode $numNodes
	} elseif {$frameType == "Sidesway_Inhibited"} {
		set controlledNode $midNode
	} else {
    	puts "frameType: $frameType not recgonized"
	}

    set dispStepSize $coarseDispStepSize

    while {$ok == 0} {

        integrator DisplacementControl $controlledNode 1 $dispStepSize
        test NormUnbalance $baseForceTolerance 10 $testOutputFlag
        algorithm Newton
        set ok [analyze 1]

        if {$ok != 0} {
            integrator DisplacementControl $controlledNode 1 [expr 0.1*$dispStepSize]
            test NormUnbalance $baseForceTolerance 10 $testOutputFlag
            algorithm Newton
            set ok [analyze 1]
        }
        if {$ok != 0} {
            integrator DisplacementControl $controlledNode 1 [expr 0.01*$dispStepSize]
            test NormUnbalance $baseForceTolerance 10 $testOutputFlag
            algorithm Newton
            set ok [analyze 1]
        }
        if {$ok != 0} {
            integrator DisplacementControl $controlledNode 1 [expr 0.1*$dispStepSize]
            test NormDispIncr $baseDisplacementTolerance 10 $testOutputFlag
            algorithm Newton
            set ok [analyze 1]
        }
        if {$ok != 0} {
            integrator DisplacementControl $controlledNode 1 [expr 0.01*$dispStepSize]
            test NormDispIncr $baseDisplacementTolerance 10 $testOutputFlag
            algorithm Newton
            set ok [analyze 1]
        }

        if {$ok == 0} {
            set lowestEigenValue [eigen -standard -symmBandLapack 1]
            puts $eigenFileId "[getTime] $lowestEigenValue"
            if { $lowestEigenValue < 0 } {
                # Analysis completed sucessfully 
                exit 1
            }
            if { [nodeDisp $controlledNode 1] > $maxDisp } {
                # Maximum deformation reached
                exit 7
            }
            incr iStep

            if { [expr $lowestEigenValue/$firstEigenValue] > 0.25 } {
                set dispStepSize $coarseDispStepSize
            } else {
                set dispStepSize $fineDispStepSize
            }

        } else {
            exit 4
        }

        if { [expr [clock seconds]-$startTime] > $maxAnalysisTime } {
            # Reached Maximum Time
            exit 8
        }   
    }
    # Step limit reached
    exit 2


# ############ DISPLACEMENT CONTROL / NON-PROPORTIONAL LOADING ############
} elseif {$analysisType == "LimitPoint_NonProportional"} {
    if {$P != 0 } { 
        # Apply gravity loads
        pattern Plain 1 Linear {
            load  $numNodes 0.0 $P 0.0  
            if { $gamma != 0 } { 
                load  101 0.0 [expr $gamma*$P] 0.0  
            } 
        }
        set ok 0
        set iStep 0
        integrator LoadControl [expr 1.0/$numStepsGravity]
        test NormUnbalance $baseForceTolerance 30 $testOutputFlag
        analysis Static
        while {$numStepsGravity > $iStep && $ok == 0} {
            set ok [analyze 1]
            if {$ok == 0} {
                set lowestEigenValue [eigen -standard -symmBandLapack 1]
                puts $eigenFileId "[getTime] $lowestEigenValue"
                incr iStep
            }   
        }
        if { $ok != 0 } {
            # Analysis failed in gravity loading
            exit 3
        }

        loadConst -time 0.0
    }

    # Apply lateral loads
    pattern Plain 2 Linear {
    	if {$frameType == "Sidesway_Uninhibited"} {
			load $numNodes 1.0 0.0 0.0  
		} elseif {$frameType == "Sidesway_Inhibited"} {
	        load 1         0.0 0.0 -1.0
        	load $numNodes 0.0 0.0 $beta
		} else {
		    puts "frameType: $frameType not recgonized"
		}
    }
    set ok 0
    set iStep 0

	if {$frameType == "Sidesway_Uninhibited"} {
		set controlledNode $numNodes
        set controlledDOF  1
	} elseif {$frameType == "Sidesway_Inhibited"} {
        set controlledNode $midNode
        set controlledDOF  1
        #set controlledNode 1
        #set controlledDOF  3
	} else {
    	puts "frameType: $frameType not recgonized"
	}

    set dispStepSize $coarseDispStepSize

    while {$ok == 0} {
        
        integrator DisplacementControl $controlledNode $controlledDOF $dispStepSize
        test NormUnbalance $baseForceTolerance 10 $testOutputFlag
        algorithm Newton
        set ok [analyze 1]
        
        if {$ok != 0} {
            integrator DisplacementControl $controlledNode $controlledDOF [expr 0.1*$dispStepSize]
            test NormUnbalance $baseForceTolerance 10 $testOutputFlag
            algorithm Newton
            set ok [analyze 1]
        }
        if {$ok != 0} {
            integrator DisplacementControl $controlledNode $controlledDOF [expr 0.01*$dispStepSize]
            test NormUnbalance $baseForceTolerance 10 $testOutputFlag
            algorithm Newton
            set ok [analyze 1]
        }
        if {$ok != 0} {
            integrator DisplacementControl $controlledNode $controlledDOF [expr 0.1*$dispStepSize]
            test NormDispIncr $baseDisplacementTolerance 10 $testOutputFlag
            algorithm Newton
            set ok [analyze 1]
        }
        if {$ok != 0} {
            integrator DisplacementControl $controlledNode $controlledDOF [expr 0.01*$dispStepSize]
            test NormDispIncr $baseDisplacementTolerance 10 $testOutputFlag
            algorithm Newton
            set ok [analyze 1]
        }
        if {$ok != 0} {
            integrator DisplacementControl $controlledNode $controlledDOF [expr 0.01*$dispStepSize]
            test NormDispIncr [expr 1e1*$baseDisplacementTolerance] 10 $testOutputFlag
            algorithm Newton  
            set ok [analyze 1]
        }
        if {$ok != 0} {
            integrator DisplacementControl $controlledNode $controlledDOF [expr 0.01*$dispStepSize]
            test NormDispIncr [expr 1e2*$baseDisplacementTolerance] 10 $testOutputFlag
            algorithm Newton  
            set ok [analyze 1]
        }
        if {$ok != 0 && $tryNumber >= 2} {
            integrator DisplacementControl $controlledNode $controlledDOF [expr 0.01*$dispStepSize]
            test NormDispIncr [expr 1e3*$baseDisplacementTolerance] 10 $testOutputFlag
            algorithm Newton  
            set ok [analyze 1]
        }
        if {$ok != 0 && $tryNumber >= 2} {
            integrator DisplacementControl $controlledNode $controlledDOF [expr 0.01*$dispStepSize]
            test NormDispIncr [expr 1e4*$baseDisplacementTolerance] 10 $testOutputFlag
            algorithm Newton  
            set ok [analyze 1]
        }
        
        if {$ok == 0} {
            set lowestEigenValue [eigen -standard -symmBandLapack 1]
            puts $eigenFileId "[getTime] $lowestEigenValue"
            if { $lowestEigenValue < 0 } {
                # Analysis completed sucessfully 
                exit 1
            }
            if { [expr abs([nodeDisp $controlledNode $controlledDOF])] > $maxDisp } {
                # Maximum deformation reached
                exit 7
            }
            incr iStep

            if { [expr $lowestEigenValue/$firstEigenValue] > 0.25 } {
                set dispStepSize $coarseDispStepSize
            } else {
                set dispStepSize $fineDispStepSize
            }

        } else {
            exit 4
        }

        if { [expr [clock seconds]-$startTime] > $maxAnalysisTime } {
            # Reached Maximum Time
            exit 8
        }   
    }
    # Step limit reached
    exit 2

# ############ LOAD CONTROL / PROPORTIONAL LOADING ############
} elseif {$analysisType == "TargetForce_Proportional"} {
    # Apply gravity loads
    pattern Plain 1 Linear {
    	if {$frameType == "Sidesway_Uninhibited"} {
            load $numNodes $H -1.0 0.0  
            if { $gamma != 0 } { 
                load  101 0.0 [expr -1*$gamma] 0.0  
            }
		} elseif {$frameType == "Sidesway_Inhibited"} {         
            load 1         0.0  0.0 [expr  -1.0*$M]
        	load $numNodes 0.0 -1.0 [expr $beta*$M]
		} else {
		    puts "frameType: $frameType not recgonized"
		}
    }

    set ok 0
    set iStep 0
    set loadStep [expr double($P)/$numStepsLateral]
    integrator LoadControl $loadStep
    test NormUnbalance $baseForceTolerance 30 $testOutputFlag
    analysis Static
    while {$numStepsLateral > $iStep && $ok == 0} {
        set ok [analyze 1]
        if {$ok != 0} {
            test NormUnbalance $baseForceTolerance 30 $testOutputFlag
            algorithm NewtonLineSearch .8
            set ok [analyze 1]
            test NormUnbalance $baseForceTolerance 30 $testOutputFlag
            algorithm Newton
        }
        if {$ok != 0} {
            test NormUnbalance $baseForceTolerance 30 $testOutputFlag
            algorithm ModifiedNewton -initial
            set ok [analyze 1]
            test NormUnbalance $baseForceTolerance 30 $testOutputFlag
            algorithm Newton
        }
        if {$ok == 0} {
            set lowestEigenValue [eigen -standard -symmBandLapack 1]
            puts $eigenFileId "[getTime] $lowestEigenValue"
            incr iStep
        }
    }

    if {$ok == 0} {
        # Analysis completed sucessfully
        exit 5
    } else {
        # Analysis failed in lateral loading
        exit 6
    }

# ############ LOAD CONTROL / NON-PROPORTIONAL LOADING ############
} elseif {$analysisType == "TargetForce_NonProportional"} {
    if { $P != 0} {
        # Apply gravity loads
        pattern Plain 1 Linear {
            load  $numNodes 0.0 $P 0.0  
            if { $gamma != 0 } { 
                load  101 0.0 [expr $gamma*$P] 0.0  
            } 
        }
        set ok 0
        set iStep 0
        integrator LoadControl [expr 1.0/$numStepsGravity]
        test NormUnbalance $baseForceTolerance 30 $testOutputFlag
        analysis Static
        while {$numStepsGravity > $iStep && $ok == 0} {
            set ok [analyze 1]
            if {$ok == 0} {
                set lowestEigenValue [eigen -standard -symmBandLapack 1]
                puts $eigenFileId "[getTime] $lowestEigenValue"
                incr iStep
            }
        }

        if { $ok != 0 } {
            # Analysis failed in gravity loading
            exit 3
        }

        loadConst -time 0.0
    }

    # Apply lateral loads
    pattern Plain 2 Linear {
    	if {$frameType == "Sidesway_Uninhibited"} {
			load $numNodes 1.0 0.0 0.0  
			set loadStep [expr double($H)/$numStepsLateral]
		} elseif {$frameType == "Sidesway_Inhibited"} {
	        load 1         0.0 0.0 -1.0
        	load $numNodes 0.0 0.0 $beta
			set loadStep [expr double($M)/$numStepsLateral]        	
		} else {
		    puts "frameType: $frameType not recgonized"
		}    
    }
    set ok 0
    set iStep 0
    integrator LoadControl $loadStep
    test NormUnbalance $baseForceTolerance 30 $testOutputFlag
    analysis Static
    while {$numStepsLateral > $iStep && $ok == 0} {
        set ok [analyze 1]
        if {$ok != 0} {
            test NormUnbalance $baseForceTolerance 30 $testOutputFlag
            algorithm NewtonLineSearch .8
            set ok [analyze 1]
            test NormUnbalance $baseForceTolerance 30 $testOutputFlag
            algorithm Newton
        }
        if {$ok != 0} {
            test NormUnbalance $baseForceTolerance 30 $testOutputFlag
            algorithm ModifiedNewton -initial
            set ok [analyze 1]
            test NormUnbalance $baseForceTolerance 30 $testOutputFlag
            algorithm Newton
        }
        if {$ok == 0} {
            set lowestEigenValue [eigen -standard -symmBandLapack 1]
            puts $eigenFileId "[getTime] $lowestEigenValue"
            incr iStep
        }
    }

    if {$ok == 0} {
        # Analysis completed sucessfully
        exit 5
    } else {
        # Analysis failed in lateral loading
        exit 6
    }

# ############ LOAD CONTROL / NON-PROPORTIONAL LOADING ############
} elseif {$analysisType == "TargetDrift_NonProportional"} {
    if { $P != 0} {
        # Apply gravity loads
        pattern Plain 1 Linear {
            load  $numNodes 0.0 $P 0.0  
            if { $gamma != 0 } { 
                load  101 0.0 [expr $gamma*$P] 0.0  
            } 
        }
        set ok 0
        set iStep 0
        integrator LoadControl [expr 1.0/$numStepsGravity]
        test NormUnbalance $baseForceTolerance 30 $testOutputFlag
        analysis Static
        while {$numStepsGravity > $iStep && $ok == 0} {
            set ok [analyze 1]
            if {$ok == 0} {
                set lowestEigenValue [eigen -standard -symmBandLapack 1]
                puts $eigenFileId "[getTime] $lowestEigenValue"
                incr iStep
            }
        }

        if { $ok != 0 } {
            # Analysis failed in gravity loading
            exit 3
        }

        loadConst -time 0.0
    }

    # Apply lateral loads
    pattern Plain 2 Linear {
    	if {$frameType == "Sidesway_Uninhibited"} {
			load $numNodes 1.0 0.0 0.0  
		} elseif {$frameType == "Sidesway_Inhibited"} {
	        puts "TargetDisplacement_NonProportional not yet implemented for Sidesway_Inhibited"      	
		} else {
		    puts "frameType: $frameType not recgonized"
		}    
    }
    set ok 0
    set iStep 0
	set controlledNode $numNodes
    set controlledDOF  1
    set dispStepSize [expr $d/double($numStepsLateral)]
    integrator DisplacementControl $controlledNode $controlledDOF $dispStepSize
    test NormUnbalance $baseForceTolerance 30 $testOutputFlag
    analysis Static
    while {$numStepsLateral > $iStep && $ok == 0} {
        set ok [analyze 1]
        if {$ok != 0} {
            test NormUnbalance $baseForceTolerance 30 $testOutputFlag
            algorithm NewtonLineSearch .8
            set ok [analyze 1]
            test NormUnbalance $baseForceTolerance 30 $testOutputFlag
            algorithm Newton
        }
        if {$ok != 0} {
            test NormUnbalance $baseForceTolerance 30 $testOutputFlag
            algorithm ModifiedNewton -initial
            set ok [analyze 1]
            test NormUnbalance $baseForceTolerance 30 $testOutputFlag
            algorithm Newton
        }
        if {$ok == 0} {
            set lowestEigenValue [eigen -standard -symmBandLapack 1]
            puts $eigenFileId "[getTime] $lowestEigenValue"
            incr iStep
        }
    }

    if {$ok == 0} {
        # Analysis completed sucessfully
        exit 5
    } else {
        # Analysis failed in lateral loading
        exit 6
    }

# ############ END OF ANALYSIS ############
} else {
    puts "analysisType: $analysisType not recgonized"
}

close $eigenFileId
