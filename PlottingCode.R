probs = matrix(0,2,2)
probs[1,1] = 0.1
probs[2,2] = 0.9
probs[1,2] = 0.3
probs[2,1] = 0.3 

sbm <- sample_sbm(n = 20, pref.matrix = probs, block.sizes = c(10,10), directed = FALSE, loops = FALSE)

converted_sbm = network::network(as_adjacency_matrix(sbm))

colour_vector = c(rep(c("tomato"), 10),rep(c("steelblue"), 10)) # these are the two clusters colour_vector[2] = "purple" # anchor point colour_vector[18] = "purple" # anchor point
colour_vector = rep(c("black"), 20)

ggnet2(converted_sbm ,color = colour_vector)

ggnet2(converted_sbm,color = colour_vector,mode="circle") # keeps layout same each time for plots 
#other options under the "Node placement" 

## Here is an idea of how to make plots from your PG output.
colours = rep("tomato",num_nodes); # default initialisation for every node
colour_options = c("tomato ","steelblue"); # this vector should be of the same length as the number of clusters you want to colour

for (i in 1:length(updated.c.bar)) {
  colours[updated.c.bar[[i]]] = colour_options[i]
}

# You can then use this colour vector in the command
ggnet2(sbm,color = colours)




