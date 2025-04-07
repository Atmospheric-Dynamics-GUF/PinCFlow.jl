function smooth!(flxwkb, state, lsmth_wkb = false)
  # do nothing
end

function smooth!(flxwkb, state, lsmth_wkb = true)

  (; sizex, sizey) = state.namelists.domain 
  (; sm_filter) = state.namelists.WKBNamelist

  if sizey == 1
    if sizex > 1 
      smooth_wkb!(flxwkb, state, sm_filter, XZ())
    else
      error("Smoothing just in z not yet implemented.")
    end
  elseif sizex == 1 
      smooth_wkb!(flxwkb, state, sm_filter, YZ())
  elseif sizex > 1
    smooth_wkb!(flxwkb, state, sm_filter, XYZ())
  end

  # TODO set boundary call still missing

end