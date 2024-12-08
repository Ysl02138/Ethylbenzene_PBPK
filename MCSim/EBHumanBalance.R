Integrate (Lsodes, 1e-9, 1e-9 , 1);
OutputFile("HumanBalance.out");

Simulation {
  
  Species =  1 ;
  BWmeas =  80;
  ExpoInduc = 1.0; 
  
  expoday = PerDose (500, 24,   0, 8);    # PerDose(exposure event, initial time, exposure duration);
  expowk  = PerDose (1.0, 168,  0, 72);   # PerDose(exposure on/off, weekly cycle (whole week/hours), initial time, exp.duration)
  expodur = PerDose (1.0, 168,  0, 384);  # PerDose(exposure on/off, study cycle (whole study period/hours), initial time, exp.duration) 
  
  PrintStep(MassBalEB,  0, 168, 1);
}
END.