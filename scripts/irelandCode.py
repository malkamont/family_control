###############
## WORKSPACE ##
###############

from scipy.stats import binomtest, skew
import graph_tool.all as gt
import matplotlib as mpl
import leidenalg as la
import igraph as ig
import pandas as pd
import numpy as np
import logging as lg
logger = lg.getLogger()
logger.setLevel("DEBUG")

###############
## FUNCTIONS ##
###############

def marginal_likelihood(g, max_neg_log = 0):
    ml = g.copy()
    if ml.is_directed():
        print("Implement directed version.")
        return
    ks = ml.strength(weights = "weight")
    td = np.sum(ks)
    for e in ml.es:
        i0, i1 = e.source, e.target
        try:
            w = e["weight"]
            ku = ks[i0]
            kv = ks[i1]
            q = td / 2.0
            p = ku * kv * 1.0 / q / q / 2.0
            p = binomtest(k = int(w), n = int(q), p = p, alternative = "greater").pvalue
            e["p_value"] = max(max_neg_log, max_neg_log if p <= 0 else p)
        except ValueError as error:
            logger.warning("warning: ValueError {}".format(str(error)))
            logger.debug("ValueError weight: {} ks[i0]: {} ks[i1]: {} td: {} p: {}".format(e["weight"], ks[i0], ks[i1], td, p))
            e["p_value"] = None
        except Exception as error:
            logger.warning("warning: Exception {}".format(str(error)))
            e["p_value"] = None
    max_sig = np.min([s for s in ml.es["p_value"] if s is not None])
    for e in ml.es:
        if e["p_value"] is None:
            e["p_value"] = max_sig
    return ml

#############
## COLOURS ##
#############

clis = ["#ff00ff", "#00ff00", "#00ffff", "#ff0000", "#ffff00", "#ff8000", "#8000ff", "#ff0080", "#80ff00", "#0080ff", "#3c4142"]
cmap = mpl.colors.ListedColormap(clis)



##############################
## ENDORSEMENT ALL ACCOUNTS ##
##############################

#graph
g = ig.Graph.Read_GraphML("/m/triton/work/malkama5/paulCoalition/irelandEndorseAllAccount.xml")
print([g.vcount(), g.ecount(), int(np.sum(g.es["weight"])), skew(g.es["weight"]), g.density()])
g = marginal_likelihood(g = g)
g.delete_edges([e.index for e in g.es if e["p_value"] >= 0.1])
g = g.connected_components().giant()
print([g.vcount(), g.ecount(), int(np.sum(g.es["weight"])), skew(g.es["weight"]), g.density()])

#cluster
op = la.Optimiser()
rs = np.round(np.arange(0.00, 1.50, 0.01), 10).tolist()
ll = [] #entropy
nc = [] #number of communities
pg = [] #graph
for r in rs:
    pz = []
    for seed in np.arange(100):
        p = la.RBConfigurationVertexPartition(g, initial_membership = None, weights = None, resolution_parameter = r)
        op.set_rng_seed(seed)
        op.optimise_partition(p, n_iterations = 10)
        pz.append(p.membership)
    gg = g.to_graph_tool(vertex_attributes = {"name": "string", "lvl": "int", "org": "string", "nam": "string", "typ": "string"})
    gt.seed_rng(1); np.random.seed(1)
    pm = gt.PartitionModeState(bs = pz, relabel = True, converge = True)
    pv = pm.get_marginal(gg)
    pb = gg.new_vp("int")
    for v in gg.vertices():
        pb[v] = np.argmax(pv[v])
    ll.append(gt.PPBlockState(gg, b = pb.a).entropy())
    nc.append(len(np.unique(pb.a)))
    gg.vp.pv = pv
    gg.vp.pb = pb
    pg.append(gg)
    print(r)

#store
g = pg[np.argmin(ll)]
df = pd.DataFrame({"cluster": g.vp.pb, "type": g.vp.typ, "level": g.vp.lvl, "user": g.vp.name, "org_id": g.vp.org, "org_name": g.vp.nam})
vc = df["cluster"].value_counts()
mp = {v:k for k, v in dict(enumerate(vc.index)).items()}
df["cluster"] = df["cluster"].map(mp)
g.vp.pb.a = df["cluster"]

df = df.sort_values(["cluster", "type", "level"])
hd = ";".join([str(e) for e in df.columns.tolist()])
np.savetxt("/m/triton/scratch/work/malkama5/paulCoalition/irelandEndorseAllAccountCluster.csv", df.values, delimiter = ";", fmt = "%s", header = hd, comments = "")
m = gt.BlockState(g, b = g.vp.pb).get_matrix().todense()

#draw
gt.seed_rng(1); np.random.seed(1)
po = gt.sfdp_layout(g, groups = g.vp.pb, gamma = 0.001, kappa = 50, max_iter = 0)
vs = g.new_vp("int")
for v in g.vertices():
    if g.vp.lvl[v] == 0:
        vs[v] = 6
    if g.vp.lvl[v] == 1:
        vs[v] = 4
    if g.vp.lvl[v] == 2:
        vs[v] = 4
    if g.vp.lvl[v] == 3:
        vs[v] = 4
gt.graph_draw(g, pos = po, vertex_size = vs, vorder = vs, vertex_fill_color = g.vp.pb, vertex_color = "#ffffff00", edge_pen_width = 0.2, vcmap = cmap, edge_color = "#d3d3d3", output = "/m/triton/scratch/work/malkama5/paulCoalition/irelandEndorseAllAccountCluster.svg")
#gt.graph_draw(g, pos = po, vertex_size = vs, vorder = vs, edge_pen_width = 0.2, vertex_color = "#00000000", vertex_shape = "pie", vertex_pie_fractions = g.vp.pv, vertex_pie_colors = clis, edge_color = "#1b1919", output = "/m/triton/scratch/work/malkama5/paulCoalition/twitterClusterAccountProbability.svg")

#directed centrality
dg = ig.Graph.Read_GraphML("/m/triton/work/malkama5/paulCoalition/irelandEndorseAllAccountDirected.xml").connected_components(mode = "weak").giant()
ideg = dg.strength(weights = "weight", mode = "in")
odeg = dg.strength(weights = "weight", mode = "out")
dg = dg.to_graph_tool(vertex_attributes = {"name": "string"}, edge_attributes = {"weight": "float"})
dg.ep.weight.a = np.log1p(dg.ep.weight.a)
pgra = gt.pagerank(dg, weight = dg.ep.weight)

d = pd.read_csv("/m/triton/scratch/work/malkama5/paulCoalition/irelandEndorseAllAccountCluster.csv", header = 0, sep = ";")
d = pd.merge(d, pd.DataFrame({"user": dg.vp.name, "ideg": ideg, "odeg": odeg, "pgra": pgra}), on = "user", how = "inner")
d["pgra"] = (d["pgra"] - d["pgra"].min()) / (d["pgra"].max() - d["pgra"].min())

clus = []
topi = []
topia = []
topo = []
topoa = []
pgra = []
pgraa = []
for c in d["cluster"].unique():
    clus.append(c)
    clus.append(c)
    clus.append(c)
    clus.append(c)
    clus.append(c)
    
    topi.append(int(d[d["cluster"] == c].sort_values("ideg", ascending = False)["ideg"].iloc[0]))
    topi.append(int(d[d["cluster"] == c].sort_values("ideg", ascending = False)["ideg"].iloc[1]))
    topi.append(int(d[d["cluster"] == c].sort_values("ideg", ascending = False)["ideg"].iloc[2]))
    topi.append(int(d[d["cluster"] == c].sort_values("ideg", ascending = False)["ideg"].iloc[3]))
    topi.append(int(d[d["cluster"] == c].sort_values("ideg", ascending = False)["ideg"].iloc[4]))
    
    topia.append(d[d["cluster"] == c].sort_values("ideg", ascending = False)["user"].iloc[0])
    topia.append(d[d["cluster"] == c].sort_values("ideg", ascending = False)["user"].iloc[1])
    topia.append(d[d["cluster"] == c].sort_values("ideg", ascending = False)["user"].iloc[2])
    topia.append(d[d["cluster"] == c].sort_values("ideg", ascending = False)["user"].iloc[3])
    topia.append(d[d["cluster"] == c].sort_values("ideg", ascending = False)["user"].iloc[4])
    
    topo.append(int(d[d["cluster"] == c].sort_values("odeg", ascending = False)["odeg"].iloc[0]))
    topo.append(int(d[d["cluster"] == c].sort_values("odeg", ascending = False)["odeg"].iloc[1]))
    topo.append(int(d[d["cluster"] == c].sort_values("odeg", ascending = False)["odeg"].iloc[2]))
    topo.append(int(d[d["cluster"] == c].sort_values("odeg", ascending = False)["odeg"].iloc[3]))
    topo.append(int(d[d["cluster"] == c].sort_values("odeg", ascending = False)["odeg"].iloc[4]))
    
    topoa.append(d[d["cluster"] == c].sort_values("odeg", ascending = False)["user"].iloc[0])
    topoa.append(d[d["cluster"] == c].sort_values("odeg", ascending = False)["user"].iloc[1])
    topoa.append(d[d["cluster"] == c].sort_values("odeg", ascending = False)["user"].iloc[2])
    topoa.append(d[d["cluster"] == c].sort_values("odeg", ascending = False)["user"].iloc[3])
    topoa.append(d[d["cluster"] == c].sort_values("odeg", ascending = False)["user"].iloc[4])
    
    pgra.append(float(d[d["cluster"] == c].sort_values("pgra", ascending = False)["pgra"].iloc[0]))
    pgra.append(float(d[d["cluster"] == c].sort_values("pgra", ascending = False)["pgra"].iloc[1]))
    pgra.append(float(d[d["cluster"] == c].sort_values("pgra", ascending = False)["pgra"].iloc[2]))
    pgra.append(float(d[d["cluster"] == c].sort_values("pgra", ascending = False)["pgra"].iloc[3]))
    pgra.append(float(d[d["cluster"] == c].sort_values("pgra", ascending = False)["pgra"].iloc[4]))
    
    pgraa.append(d[d["cluster"] == c].sort_values("pgra", ascending = False)["user"].iloc[0])
    pgraa.append(d[d["cluster"] == c].sort_values("pgra", ascending = False)["user"].iloc[1])
    pgraa.append(d[d["cluster"] == c].sort_values("pgra", ascending = False)["user"].iloc[2])
    pgraa.append(d[d["cluster"] == c].sort_values("pgra", ascending = False)["user"].iloc[3])
    pgraa.append(d[d["cluster"] == c].sort_values("pgra", ascending = False)["user"].iloc[4])
dcn = pd.DataFrame({"cluster": clus, "top_retweeted_account": topoa, "top_retweeted_count": topo, "top_retweeter_account": topia, "top_retweeter_count": topi, "top_pagerank_account": pgraa, "top_pagerank_value": pgra}).round(2)
dcn.sort_values(["cluster", "top_pagerank_value"])

#edge bundle
ebun = pd.merge(pd.DataFrame({"user": g.vp.name}), d, on = "user", how = "left")
pgra = g.new_vp("float")
for i, v in enumerate(g.vertices()):
    pgra[v] = ebun["pgra"].iloc[i]

ebundle = gt.NestedBlockState(g, bs = [g.vp.pb, g.vp.pb])
gt.draw_hierarchy(ebundle, beta = 1.5, vertex_color = "#00000000", vcmap = cmap, rel_order = pgra, vorder = pgra, edge_pen_width = 0.1, vertex_size = gt.prop_to_size(pgra, 2, 2), hvertex_size = 0, hedge_pen_width = 0, hedge_marker_size = 0, vertex_text_color = "#000000", vertex_text = g.vp.name, vertex_text_position = "centered", vertex_text_offset = [0.01, 0], output = "/m/triton/scratch/work/malkama5/paulCoalition/irelandEndorseAllAccountClusterEdgeBundle.svg")



#################
## CO-ADVOCACY ##
#################

#graph
g = ig.Graph.Read_GraphML("/m/triton/work/malkama5/paulCoalition/irelandCoAdvocacy.xml")
print([g.vcount(), g.ecount(), int(np.sum(g.es["weight"])), skew(g.es["weight"]), g.density()])
g = marginal_likelihood(g = g)
g.delete_edges([e.index for e in g.es if e["p_value"] >= 0.1])
g = g.connected_components().giant()
print([g.vcount(), g.ecount(), int(np.sum(g.es["weight"])), skew(g.es["weight"]), g.density()])

#cluster
op = la.Optimiser()
rs = np.round(np.arange(0.00, 1.50, 0.01), 10).tolist()
ll = [] #entropy
nc = [] #number of communities
pg = [] #graph
for r in rs:
    pz = []
    for seed in np.arange(100):
        p = la.RBConfigurationVertexPartition(g, initial_membership = None, weights = None, resolution_parameter = r)
        op.set_rng_seed(seed)
        op.optimise_partition(p, n_iterations = 10)
        pz.append(p.membership)
    gg = g.to_graph_tool(vertex_attributes = {"name": "string", "nam": "string", "typ": "string"})
    gt.seed_rng(1); np.random.seed(1)
    pm = gt.PartitionModeState(bs = pz, relabel = True, converge = True)
    pv = pm.get_marginal(gg)
    pb = gg.new_vp("int")
    for v in gg.vertices():
        pb[v] = np.argmax(pv[v])
    ll.append(gt.PPBlockState(gg, b = pb.a).entropy())
    nc.append(len(np.unique(pb.a)))
    gg.vp.pv = pv
    gg.vp.pb = pb
    pg.append(gg)
    print(r)

#store
g = pg[np.argmin(ll)]
df = pd.DataFrame({"cluster": g.vp.pb, "type": g.vp.typ, "org_id": g.vp.name, "org_name": g.vp.nam})
vc = df["cluster"].value_counts()
mp = {v:k for k, v in dict(enumerate(vc.index)).items()}
df["cluster"] = df["cluster"].map(mp)
g.vp.pb.a = df["cluster"]

df = df.sort_values(["cluster", "type"])
hd = ";".join([str(e) for e in df.columns.tolist()])
np.savetxt("/m/triton/scratch/work/malkama5/paulCoalition/irelandCoAdvocacyCluster.csv", df.values, delimiter = ";", fmt = "%s", header = hd, comments = "")
m = gt.BlockState(g, b = g.vp.pb).get_matrix().todense()
dsur = df.copy()

#coordinates
gt.seed_rng(1); np.random.seed(1)
po = gt.sfdp_layout(g, groups = g.vp.pb, gamma = 0.00001, kappa = 20)
#s.draw(pos = p, vertex_text = g.vp.name)

#draw
#pp = pd.merge(pd.DataFrame({"vertex": g.vp.name}), P, on = "vertex", how = "left")
#po = g.new_vp("vector<double>")
#for i, v in enumerate(g.vertices()):
#    po[v] = pp["coord"][i]
gt.graph_draw(g, pos = po, vertex_size = 18, vertex_fill_color = g.vp.pb, vertex_color = "#ffffff00", edge_pen_width = 0.3, vcmap = cmap, edge_color = "#d3d3d3", vertex_text = g.vp.nam, vertex_font_size = 0, vertex_text_out_color = "#ffffff", vertex_text_color = "#1b1919", vertex_text_rotation = 6.1, vertex_text_position = 0, vertex_text_offset = [0.001, 0.000], vertex_text_out_width = 0.00, output = "/m/triton/scratch/work/malkama5/paulCoalition/irelandCoAdvocacyCluster.svg")
#gt.graph_draw(g, pos = po, vertex_size = 10, edge_pen_width = 0.3, vertex_color = "#00000000", vertex_shape = "pie", vertex_pie_fractions = g.vp.pv, vertex_pie_colors = clis, edge_color = "#1b1919", output = "/m/triton/scratch/work/malkama5/paulCoalition/surveyClusterOrganisationProbability.svg")

#complete centrality
cg = ig.Graph.Read_GraphML("/m/triton/work/malkama5/paulCoalition/irelandCoAdvocacy.xml")
cg.es["weight"] = np.log1p(cg.es["weight"])
cg = cg.to_graph_tool(vertex_attributes = {"name": "string", "nam": "string", "typ": "string"}, edge_attributes = {"weight": "float"})
ev = gt.eigenvector(cg, weight = cg.ep.weight)[1]
cl = gt.closeness(cg, weight = cg.ep.weight)

d = pd.DataFrame({"cluster": g.vp.pb, "org": g.vp.name, "eigv": ev, "clos": cl})
d["eigv"] = (d["eigv"] - d["eigv"].min()) / (d["eigv"].max() - d["eigv"].min())
d["clos"] = (d["clos"] - d["clos"].min()) / (d["clos"].max() - d["clos"].min())
pd.concat([
d[d["cluster"] == 0].sort_values("eigv", ascending = False).head(3),
d[d["cluster"] == 1].sort_values("eigv", ascending = False).head(3)], ignore_index = True, axis = 0)

EV = g.new_vp("float")
for i, v in enumerate(g.vertices()):
    EV[v] = ev[i]

#edge bundle
ebundle = gt.NestedBlockState(g, bs = [g.vp.pb, g.vp.pb])
gt.draw_hierarchy(ebundle, beta = 0.8, vertex_color = "#00000000", vcmap = cmap, rel_order = EV, vorder = EV, edge_pen_width = 0.4, vertex_size = gt.prop_to_size(EV, 8, 8), hvertex_size = 0, hedge_pen_width = 0, hedge_marker_size = 0, vertex_text_color = "#000000", vertex_text = g.vp.nam, vertex_text_position = "centered", vertex_text_offset = [0.02, 0], output = "/m/triton/scratch/work/malkama5/paulCoalition/irelandCoAdvocacyClusterEdgeBundle.svg")



#################################################
## ENDORSEMENT CLUSTERS IN CO-ADVOCACY NETWORK ##
#################################################

#assignments
ia = pd.read_csv("/m/triton/scratch/work/malkama5/paulCoalition/irelandAttributes.csv", delimiter = ";", header = 0)
aa = pd.read_csv("/m/triton/scratch/work/malkama5/paulCoalition/twitterAccountAssignment.csv", delimiter = ";", header = 0)
df = pd.merge(pd.DataFrame({"org": g.vp.name}), ia, on = "org", how = "left")
df = pd.merge(df, aa, on = "name", how = "left")
df.iloc[:, 3:13] = df.iloc[:, 3:13].fillna(0)
df["author_majority"] = df["author_majority"].fillna(999)

pv = g.new_vp("vector<float>")
pm = g.new_vp("vector<float>")
pb = g.new_vp("int")
for v in np.arange(g.num_vertices()):
    x = df.iloc[v, 3:13].tolist()
    if np.sum(x) > 0:
        x = x / np.sum(x) * 100
        pv[v] = x
    else:
        x = x + [100]
        pv[v] = x
        x = x[:len(x) - 1]
    m = int(df.loc[v, "author_majority"])
    x = [0] * len(x)
    if m != 999:
        x[m] = 100
    else:
        x = x + [100]
    pm[v] = x
    pb[v] = np.argmax(x)

gt.graph_draw(g, pos = po, vertex_size = 18, edge_pen_width = 0.3, vertex_color = "#ffffff00", vertex_shape = "pie", vertex_pie_fractions = pm, vertex_pie_colors = clis, edge_color = "#d3d3d3", output = "/m/triton/scratch/work/malkama5/paulCoalition/irelandEndorseCoAdvocacyMajority.svg")

gt.graph_draw(g, pos = po, vertex_size = 18, edge_pen_width = 0.3, vertex_color = "#ffffff00", vertex_shape = "pie", vertex_pie_fractions = pv, vertex_pie_colors = clis, edge_color = "#d3d3d3", output = "/m/triton/scratch/work/malkama5/paulCoalition/irelandEndorseCoAdvocacyMarginal.svg")

df = pd.DataFrame({"org_id": df["org"], "cluster_author_majority": pb, "cluster_marginal": pv})
df = df[df["cluster_author_majority"] < df["cluster_author_majority"].max()]
hd = ";".join([str(e) for e in df.columns.tolist()])
np.savetxt("/m/triton/scratch/work/malkama5/paulCoalition/irelandEndorseCoAdvocacy.csv", df.values, delimiter = ";", fmt = "%s", header = hd, comments = "")
m = gt.BlockState(g, b = pb).get_matrix().todense()
ddis = df.copy()

########################
## MUTUAL INFORMATION ##
########################

mi = pd.merge(dsur.loc[:, ["org_id", "cluster"]], ddis.loc[:, ["org_id", "cluster_author_majority"]], on = "org_id", how = "inner")
mi = gt.reduced_mutual_information(mi["cluster"], mi["cluster_author_majority"], norm = True)
print(mi)

####################################
## ENDORSEMENT AGGREGATE ACCOUNTS ##
####################################

#graph
g = ig.Graph.Read_GraphML("/m/triton/work/malkama5/paulCoalition/irelandEndorseAggregate.xml").connected_components().giant()
print([g.vcount(), g.ecount(), int(np.sum(g.es["weight"])), skew(g.es["weight"]), g.density()])
g = marginal_likelihood(g = g)
g.delete_edges([e.index for e in g.es if e["p_value"] >= 0.1])
g = g.connected_components().giant()
print([g.vcount(), g.ecount(), int(np.sum(g.es["weight"])), skew(g.es["weight"]), g.density()])

#cluster
op = la.Optimiser()
rs = np.round(np.arange(0.00, 1.50, 0.01), 10).tolist()
ll = [] #entropy
nc = [] #number of communities
pg = [] #graph
for r in rs:
    pz = []
    for seed in np.arange(100):
        p = la.RBConfigurationVertexPartition(g, initial_membership = None, weights = None, resolution_parameter = r)
        op.set_rng_seed(seed)
        op.optimise_partition(p, n_iterations = 10)
        pz.append(p.membership)
    gg = g.to_graph_tool(vertex_attributes = {"name": "string", "nam": "string", "typ": "string"})
    gt.seed_rng(1); np.random.seed(1)
    pm = gt.PartitionModeState(bs = pz, relabel = True, converge = True)
    pv = pm.get_marginal(gg)
    pb = gg.new_vp("int")
    for v in gg.vertices():
        pb[v] = np.argmax(pv[v])
    ll.append(gt.PPBlockState(gg, b = pb.a).entropy())
    nc.append(len(np.unique(pb.a)))
    gg.vp.pv = pv
    gg.vp.pb = pb
    pg.append(gg)
    print(r)

#store
g = pg[np.argmin(ll)]
df = pd.DataFrame({"cluster": g.vp.pb, "type": g.vp.typ, "org_id": g.vp.name, "org_name": g.vp.nam})
vc = df["cluster"].value_counts()
mp = {v:k for k, v in dict(enumerate(vc.index)).items()}
df["cluster"] = df["cluster"].map(mp)
g.vp.pb.a = df["cluster"]

df = df.sort_values(["cluster", "type"])
hd = ";".join([str(e) for e in df.columns.tolist()])
np.savetxt("/m/triton/scratch/work/malkama5/paulCoalition/irelandEndorseAggregateCluster.csv", df.values, delimiter = ";", fmt = "%s", header = hd, comments = "")
m = gt.BlockState(g, b = g.vp.pb).get_matrix().todense()
dagg = df.copy()

#rmi
mi = pd.merge(dsur.loc[:, ["org_id", "cluster"]], dagg.loc[:, ["org_id", "cluster"]], on = "org_id", how = "inner")
mi = gt.reduced_mutual_information(mi["cluster_x"], mi["cluster_y"], norm = True)
print(mi)
