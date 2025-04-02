function meanflow_interpolation(
    flwlbd,
    flwlbu,
    flwlfd,
    flwlfu,
    zlbd,
    zlbu,
    zlfd,
    zlfu,
    xl,
    xr,
    xlc,
    yf,
    yb,
    ylc,
    zlc,
  )
    if sizex == 1
      flwbd = flwlbd
      flwbu = flwlbu
  
      flwfd = flwlfd
      flwfu = flwlfu
  
      zbd = zlbd
      zbu = zlbu
  
      zfd = zlfd
      zfu = zlfu
    else
      if xr < xl
        error("ERROR IN MEANFLOW: xr =", xr, "< xl =", xl)
      elseif xr == xl
        factor = 0.0
      elseif xlc > xr
        factor = 0.0
      elseif xlc > xl
        factor = (xr - xlc) / dx
      else
        factor = 1.0
      end
  
      flwbd = factor * flwlbd + (1.0 - factor) * flwrbd
      flwbu = factor * flwlbu + (1.0 - factor) * flwrbu
  
      flwfd = factor * flwlfd + (1.0 - factor) * flwrfd
      flwfu = factor * flwlfu + (1.0 - factor) * flwrfu
  
      zbd = factor * zlbd + (1.0 - factor) * zrbd
      zbu = factor * zlbu + (1.0 - factor) * zrbu
  
      zfd = factor * zlfd + (1.0 - factor) * zrfd
      zfu = factor * zlfu + (1.0 - factor) * zrfu
    end
  
    # interpolation in y
  
    if sizey == 1
      flwd = flwbd
      flwu = flwbu
  
      zd = zbd
      zu = zbu
    else
      if yf < yb
        error("ERROR IN MEANFLOW: yf =", yf, "< yb =", yb)
      elseif yf == yb
        factor = 0.0
      elseif ylc > yf
        factor = 0.0
      elseif ylc > yb
        factor = (yf - ylc) / dy
      else
        factor = 1.0
      end
  
      flwd = factor * flwbd + (1.0 - factor) * flwfd
      flwu = factor * flwbu + (1.0 - factor) * flwfu
  
      zd = factor * zbd + (1.0 - factor) * zfd
      zu = factor * zbu + (1.0 - factor) * zfu
    end
  
    # interpolation in z
    if zu < zd
      error("ERROR IN MEANFLOW: zu =", zu, "< zd =", zd)
    elseif zu == zd
      factor = 0.0
    elseif zlc > zu
      factor = 0.0
    elseif zlc > zd
      factor = (zu - zlc) / (zu - zd)
    else
      factor = 1.0
    end
  
    flw = factor * flwd + (1.0 - factor) * flwu
    return flw
  end