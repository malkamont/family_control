###############
## WORKSPACE ##
###############

library(igraph)

############
## X DATA ##
############

d = readRDS("/m/triton/work/malkama5/paulCoalition/irelandData.rds")
d = subset(d, !is.na(ref_org))
d = subset(d, type == "retweeted" & created_at >= as.Date("2017-12-16") & created_at <= as.Date("2021-12-15")) #four years
f = grep("biodiversity|carbon|clean energy|climate|coal|energy|fossil fuel|fracking|fridaysforfuture|global warming|greenhouse gas|ghg|gas|heatwave|mass extinction|methane|net zero|oil|paris accord|paris agreement|renewable energy|sea level|solar|sustainability|warming|wind energy|wind power", d$text, ignore.case = TRUE)
d = d[f, 1:ncol(d)]
d$author_sn_roster = paste0(d$org, "_", d$author_sn_roster)
d$ref_author_sn_roster = paste0(d$ref_org, "_", d$ref_author_sn_roster)

#######################################
## ENDORSEMENT ALL ACCOUNTS DIRECTED ##
#######################################

#graph
g = graph_from_data_frame(d[c(which(names(d) == "ref_author_sn_roster"), which(names(d) == "author_sn_roster"))], directed = TRUE)
E(g)$weight = 1
g = simplify(g, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = "sum")

#attributes
a = read.table("/m/triton/scratch/work/malkama5/paulCoalition/irelandAccounts.csv", header = TRUE, sep = ",", colClasses = "character")
m = merge(x = data.frame("actor" = V(g)$name), y = data.frame("actor" = paste0(a$org, "_", a$sn), "org" = a$org, "level" = a$level), by = "actor", all.x = TRUE, sort = FALSE)
m$index = 1:nrow(m)

a = read.table("/m/triton/scratch/work/malkama5/paulCoalition/irelandAttributes.csv", header = TRUE, sep = ";", colClasses = "character")
m = merge(x = m, y = a, by = "org", all.x = TRUE, sort = FALSE)
m = m[order(m$index), 1:ncol(m)]

V(g)$org = m$org
V(g)$lvl = m$level
V(g)$typ = m$type
V(g)$nam = m$name

#store
write_graph(g, "/m/triton/scratch/work/malkama5/paulCoalition/irelandEndorseAllAccountDirected.xml", format = "graphml")

##############################
## ENDORSEMENT ALL ACCOUNTS ##
##############################

#graph
g = graph_from_data_frame(d[c(which(names(d) == "ref_author_sn_roster"), which(names(d) == "author_sn_roster"))], directed = FALSE)
E(g)$weight = 1
g = simplify(g, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = "sum")
g = delete_vertices(g, V(g)[strength(g) < 2]) #activity check

#attributes
a = read.table("/m/triton/scratch/work/malkama5/paulCoalition/irelandAccounts.csv", header = TRUE, sep = ",", colClasses = "character")
m = merge(x = data.frame("actor" = V(g)$name), y = data.frame("actor" = paste0(a$org, "_", a$sn), "org" = a$org, "level" = a$level), by = "actor", all.x = TRUE, sort = FALSE)
m$index = 1:nrow(m)

a = read.table("/m/triton/scratch/work/malkama5/paulCoalition/irelandAttributes.csv", header = TRUE, sep = ";", colClasses = "character")
m = merge(x = m, y = a, by = "org", all.x = TRUE, sort = FALSE)
m = m[order(m$index), 1:ncol(m)]

V(g)$org = m$org
V(g)$lvl = m$level
V(g)$typ = m$type
V(g)$nam = m$name

#store
write_graph(g, "/m/triton/scratch/work/malkama5/paulCoalition/irelandEndorseAllAccount.xml", format = "graphml")

###########################
## ENDORSEMENT AGGREGATE ##
###########################

#graph
g = graph_from_data_frame(d[c(which(names(d) == "ref_org"), which(names(d) == "org"))], directed = FALSE) #aggregate by organisation
E(g)$weight = 1
g = simplify(g, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = "sum")
g = delete_vertices(g, V(g)[strength(g) < 2]) #activity check

#attributes
a = read.table("/m/triton/scratch/work/malkama5/paulCoalition/irelandAccounts.csv", header = TRUE, sep = ",", colClasses = "character")
m = merge(x = data.frame("actor" = V(g)$name), y = data.frame("actor" = paste0(a$org, "_", a$sn), "org" = a$org, "level" = a$level), by = "actor", all.x = TRUE, sort = FALSE)
m$index = 1:nrow(m)

a = read.table("/m/triton/scratch/work/malkama5/paulCoalition/irelandAttributes.csv", header = TRUE, sep = ";", colClasses = "character")
m = merge(x = m, y = a, by = "org", all.x = TRUE, sort = FALSE)
m = m[order(m$index), 1:ncol(m)]

V(g)$org = m$org
V(g)$lvl = m$level
V(g)$typ = m$type
V(g)$nam = m$name

#store
write_graph(g, "/m/triton/scratch/work/malkama5/paulCoalition/irelandEndorseAggregate.xml", format = "graphml")

##########################################
## ENDORSEMENT COLLECTIVE MAIN ACCOUNTS ##
##########################################

#graph
g = graph_from_data_frame(d[c(which(names(d) == "ref_author_sn_roster"), which(names(d) == "author_sn_roster"))], directed = FALSE)
E(g)$weight = 1
g = simplify(g, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = "sum")
g = delete_vertices(g, V(g)[strength(g) < 2]) #activity check

#attributes
a = read.table("/m/triton/scratch/work/malkama5/paulCoalition/irelandAccounts.csv", header = TRUE, sep = ",", colClasses = "character")
m = merge(x = data.frame("actor" = V(g)$name), y = data.frame("actor" = paste0(a$org, "_", a$sn), "org" = a$org, "level" = a$level), by = "actor", all.x = TRUE, sort = FALSE)
m$index = 1:nrow(m)

a = read.table("/m/triton/scratch/work/malkama5/paulCoalition/irelandAttributes.csv", header = TRUE, sep = ";", colClasses = "character")
m = merge(x = m, y = a, by = "org", all.x = TRUE, sort = FALSE)
m = m[order(m$index), 1:ncol(m)]

V(g)$org = m$org
V(g)$lvl = m$level
V(g)$typ = m$type
V(g)$nam = m$name

#purge
g = delete_vertices(g, V(g)[V(g)$lvl != 0])

#store
write_graph(g, "/m/triton/scratch/work/malkama5/paulCoalition/irelandEndorseCollectiveMain.xml", format = "graphml")

#################
## SURVEY DATA ##
#################

be = as.matrix(read.table("/m/triton/scratch/work/malkama5/paulCoalition/irelandBeliefs.csv", header = TRUE, sep = ";", row.names = 1))
de = psych::describe(be)
de = de[order(de$sd, decreasing = TRUE), 1:ncol(de)]
be = be[1:nrow(be), de$vars[de$sd > 1]]
be[be < 3] = 1 #nay
be[be == 3] = 0
be[be > 3] = 2 #yea
bb = be
be[be != 1] = 0
bb[bb != 2] = 0
be = cbind(be, bb)
be[be == 2] = 1
rowSums(be)[rowSums(be) == 0]
be = be %*% t(be)
diag(be) = 0

co = as.matrix(read.table("/m/triton/scratch/work/malkama5/paulCoalition/irelandCollaboration.csv", header = TRUE, sep = ";", row.names = 1))
co = graph_from_adjacency_matrix(co, mode = "directed", weighted = NULL, diag = FALSE)
co = as.undirected(co, mode = "each")
E(co)$weight = 1
co = simplify(co, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = "sum")
co = as_adjacency_matrix(co, attr = "weight", sparse = FALSE)
co = co[rownames(be), rownames(be)]
co[co > 1] = 1 #disregard reciprocity

for (row in 1:nrow(co)){
    for (col in 1:ncol(co)){
        if (co[row, col] == 1){
            co[row, col] = be[row, col]}}} 

#################
## CO-ADVOCACY ##
#################

#graph
g = graph_from_adjacency_matrix(co, mode = "upper", weighted = TRUE, diag = FALSE)

#attributes
m = data.frame("org" = V(g)$name)
m$index = 1:nrow(m)

a = read.table("/m/triton/scratch/work/malkama5/paulCoalition/irelandAttributes.csv", header = TRUE, sep = ";", colClasses = "character")
m = merge(x = m, y = a, by = "org", all.x = TRUE, sort = FALSE)
m = m[order(m$index), 1:ncol(m)]

V(g)$typ = m$type
V(g)$nam = m$name

#store
write_graph(g, "/m/triton/scratch/work/malkama5/paulCoalition/irelandCoAdvocacy.xml", format = "graphml")
