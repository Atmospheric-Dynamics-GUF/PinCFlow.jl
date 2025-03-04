using PinCFlow_dev

semi = initialize_values(5, 1, 5, 3, 3, 3, 0, 60000, 0, 40000, 0, 20000)
initialize_atmosphere!(semi)
dt = 5.0

compute_fluxes!(semi);

#reconstruct_rhop!(semi);
