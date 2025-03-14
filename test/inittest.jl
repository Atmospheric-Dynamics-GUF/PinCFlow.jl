using PinCFlow

semi = initialize_values(150, 1, 50, 3, 3, 3, 0, 60000, 0, 40000, 0, 20000)

dt = 30.0

pincflow(semi, dt);
