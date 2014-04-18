print_formatted(fmt, args...) = @eval @printf($fmt, $(args...))

sprint_formatted(fmt, args...) = @eval @sprintf($fmt, $(args...))

function print_array(x :: Vector; border=3, fmt :: String = "%8.1e")
  n = length(x)
  b1, b2 = (n > 2 * border) ? (border, n-border+1) : (n, n+1)
  s = "["
  for k = 1 : b1
    s *= sprint_formatted(fmt, x[k])
    s *= " "
  end
  if n > 2 * border
    s *= "..."
  end
  for k = b2 : n
    s *= " "
    s *= sprint_formatted(fmt, x[k])
  end
  s *= "]"
  println(s)
end
