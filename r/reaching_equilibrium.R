target = abs(eigen(t(matrix(all_egs[sr == 'Lockdown 2']$cms, nrow=9)))$vectors[1,])

initial = abs(eigen(t(matrix(all_egs[sr == 'Lockdown 3']$cms, nrow=9)))$vectors[1,])

transmat = t(matrix(all_egs[sr == 'Lockdown 2']$cms, nrow=9))


active = initial 

for(i in 1:5){
  
  active =  active %*% transmat
  
  print(active / sum(active) )
  
}


target / sum(target)
initial / sum(initial)
