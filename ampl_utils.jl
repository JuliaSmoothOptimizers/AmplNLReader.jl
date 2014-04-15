function print_array(x :: Array{Float64,1}; border=3, fmt="%8.1e")
  @eval print_formatted(x) = @sprintf($fmt, x)
  n = length(x)
  b1, b2 = (n > 2 * border) ? (border, n-border+1) : (n, n+1)
  s = "["
  for k = 1 : b1
    s *= print_formatted(x[k])
    s *= " "
  end
  if n > 2 * border
    s *= "..."
  end
  for k = b2 : n
    s *= " "
    s *= print_formatted(x[k])
  end
  s *= "]"
  println(s)
end
