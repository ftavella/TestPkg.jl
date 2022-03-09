module TestPkg

using DifferentialEquations

export oscillator, solve_ode

function oscillator(u,p,t)
  x, y, z = u
  A, B, G, K, N = p
  du1 = B[1]*(A[1] - x) + G[1,3] * (1.0 - x) * abs(z)^N[1,3] / (abs(z)^N[1,3] + K[1,3]^N[1,3])
  du2 = B[2]*(A[2] - y) + G[2,1] * (1.0 - y) * abs(x)^N[2,1] / (abs(x)^N[2,1] + K[2,1]^N[2,1])
  du3 = B[3]*(A[3] - z) - G[3,2] * z * abs(y)^N[3,2] / (abs(y)^N[3,2] + K[3,2]^N[3,2])
  [du1,du2,du3]
end

function solve_ode(prob)
  sol = solve(prob, isoutofdomain=(y,p,t)->any(x->x<0,y))
  return sol
end

end
