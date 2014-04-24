print_formatted(fmt, args...) = @eval @printf($fmt, $(args...))

sprint_formatted(fmt, args...) = @eval @sprintf($fmt, $(args...))
