# Community detection for stochastic block model with unbalanced cluster sizes

SDP formulation:
max <A,X>
s.t. X>=0, 0<=X_ij<=1/m_min, X1=1, tr(X)=k

where A is the adjacency matrix and k is number of clusters, m_min is the minimal cluster size.

# Reference
Yan, Bowei, and Purnamrita Sarkar. "Convex Relaxation for Community Detection with Covariates." arXiv preprint arXiv:1607.02675 (2016).

Yan, Bowei, Purnamrita Sarkar, and Xiuyuan Cheng. "Exact Recovery of Number of Blocks in Blockmodels." arXiv preprint arXiv:1705.08580 (2017).
