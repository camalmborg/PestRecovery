# array/parallel jobs test script
vec <- seq(from = 1,
           to = 20,
           by = 2)

for (i in 1:length(vec)){
  print(vec[i] + 1)
}
