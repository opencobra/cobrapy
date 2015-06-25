
# # Gapfillling
# 
# GrowMatch and SMILEY are gap-filling algorithms, which try to to make the
# minimal number of changes to a model and allow it to simulate growth. For
# more information, see [Kumar et
# al.](http://dx.doi.org/10.1371/journal.pcbi.1000308). Please note that these
# algorithms are Mixed-Integer Linear Programs, which need solvers such as
# gurobi or cplex to function correctly.



# In this model D-Fructose-6-phosphate is an essential metabolite. We will
# remove all the reactions using it, and at them to a separate model.



# Now, because of these gaps, the model won't grow.



# We will use GrowMatch to add back the minimal number of reactions from this
# set of "universal" reactions (in this case just the ones we removed) to allow
# it to grow.



# We can obtain multiple possible reaction sets by having the algorithm go
# through multiple iterations.


# Prints:
# ---- Run 1 ----
# GF6PTA
# TKT2_reverse
# MAN6PI_reverse
# PGI_reverse
# F6PA_reverse
# ---- Run 2 ----
# TALA
# F6PP
# FBP
# GF6PTA
# MAN6PI_reverse
# ---- Run 3 ----
# GF6PTA
# TKT2_reverse
# MAN6PI_reverse
# PGI_reverse
# F6PA_reverse
# ---- Run 4 ----
# TALA
# F6PP
# FBP
# GF6PTA
# MAN6PI_reverse
