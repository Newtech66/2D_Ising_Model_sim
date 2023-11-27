echo off
set arg1=%1
shift
call "C:\Users\Mainak\OneDrive\Documents\Academic\Year 4\Statistical Physics\2D_Ising_Model_sim\CPP\magnetisation.exe" %arg1%
call "C:\Users\Mainak\OneDrive\Documents\Academic\Year 4\Statistical Physics\2D_Ising_Model_sim\CPP\susceptibility.exe" %arg1%
call "C:\Users\Mainak\OneDrive\Documents\Academic\Year 4\Statistical Physics\2D_Ising_Model_sim\CPP\specific_heat.exe" %arg1%
