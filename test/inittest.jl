using PinCFlow_dev

semi = initialize_values(300, 1, 100, 3, 3, 3, 0, 60000, 0, 40000, 0, 20000)

dt = 1.0

pincflow(semi, dt)