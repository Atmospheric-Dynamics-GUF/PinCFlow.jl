function set_integrals_to_zero!(integrals::Integrals)
  for field in fieldnames(Integrals)
    getfield(integrals, field) .= 0.0
  end

  return
end