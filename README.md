# BTTree-Reduce

This is a Magma package for computations relating to the Bruhat-Tits tree over a finite extension of a p-adic field.

The functions will be described in more detail in my upcoming thesis.

## Main intrinsics

### Creation of trees and vertices

`BruhatTitsTree(field::FldPad) -> BTTree`

`BruhatTitsTree(ring::RngPad) -> BTTree`

Create a Bruhat-Tits tree from a p-adic ring or field.

`BTTVertex(tree::BTTree, expansion::FldPadElt, precision::RngIntElt) -> BTTVert`

`BTTVertex(tree::BTTree, expansion::RngPadElt, precision::RngIntElt) -> BTTVert`

`BTTVertex(tree::BTTree, expansion::RngIntElt, precision::RngIntElt) -> BTTVert`

Create a vertex as a p-adic approximation.

`BTTVertex(tree::BTTree, matrix::AlgMatElt[FldPad]) -> BTTVert`

Convert a matrix to a Bruhat-Tits tree vertex.

### Properties of trees and vertices

`intrinsic Field(tree::BTTree) -> FldPad`

Get the underlying field of the tree.

`Prime(tree::BTTree) -> RngIntElt`

`Uniformizer(tree::BTTree) -> RngIntElt`

Get a uniformizer of the tree.

`Origin(tree::BTTree) -> BTTVert`

Get the 'origin' of the tree, that is the vertex corresponding to the identity matrix

`Degree(tree:BTTree) -> RngIntElt`

Return the degree of the tree.

`Matrix(v::BTTVert) -> AlgMatElt[FldPad]`

Convert a vertex to a matrix in vertex normal form.

`Field(vertex::BTTVert) -> FldPad`

Return the underlying field of a vertex.

`Prime(vertex::BTTVert) -> RngIntElt`

Return the prime used for the valuation.

`Uniformizer(vertex::BTTVert) -> RngIntElt`

Return the prime used for the valuation.

`Parent(vertex::BTTVert) -> BTTree`

Return the Bruhat-Tits tree containing the vertex.

`Expansion(v::BTTVert) -> FldPadElt`

Return the expansion component of v.

`Precision(v::BTTVert) -> RngIntElt`

Return the precision component of v.

### Operations on vertices

`Distance(v::BTTVert, w::BTTVert) -> RngIntElt`

Return the distance between two vertices.

`Neighbor(v::BTTVert, x::ModTupFldElt[FldFin]) -> BTTVert`

`Neighbor(v::BTTVert, x::RngElt) -> BTTVert`

Return the neighbor w of v such that [1, x] is in the subspace v/w of k^2, where k is the residue field.

`Neighbor(v::BTTVert, x::Infty) -> BTTVert`

Return the neighbor w of v such that [0, 1] is in the subspace v/w of k^2, where k is the residue field

`TypeOfNeighbor(v::BTTVert, w::BTTVert) -> ModTupFldElt[FldFin]`

Given neighbouring vertices v and w, give v/w (as a rank 2 vector over the residue field) in Echelon form
  (ie. [1, u] or [0, 1], where u is in the residue field).

Return the neighbours of v.

`IsometryBetweenAxes(tree::BTTree, v::ModTupFldElt[FldPad], w::ModTupFldElt[FldPad]) -> AlgMatElt[FldPad]`

Return an isometry of translation length 1 inducing a path from v to w.

`IsometryBetweenVertices(v::BTTVert, w::BTTVert) -> SeqEnum[BTTVert]`

Return an isometry of translation length 1 inducing a path from v to w.

`Path(v::BTTVert, w::BTTVert) -> SeqEnum[BTTVert]`

Return the path between two vertices.

`IsInMinTranslationSet(A::AlgMatElt[FldPad], v::BTTVert) -> BoolElt`

`IsInMinTranslationSetBoundary(A::AlgMatElt[FldPad], v::ModTupFldElt[FldPad]) -> BoolElt`

Return true if v is in the minimum translation set of A.

### Properties of actions on the Bruhat-Tits tree

`'*'(v::BTTVert, mat::AlgMatElt) -> BTTVert`

Apply the action of a matrix to a vertex.

`μ(A::AlgMatElt) -> FldRatElt`

Return μ(A). Roughly speaking this is how far the characteristic polynomial of A is from being 'reduced'.

`TranslationLength(A::AlgMatElt) -> RngIntElt`

Return the translation length of the action.

`IsElliptic(mat::AlgMatElt) -> BoolElt`

Return true if the isometry determined by mat is elliptic.

`IsHyperbolic(mat::AlgMatElt) -> BoolElt`

Return true if the isometry determined by mat is hyperbolic.

`TranslationAxisBoundary(A::AlgMatElt[FldPad]) -> Tup`

Return the points on the boundary of the translation axis of the hyperbolic isometry determined by A.

`Intersects(A::AlgMatElt[FldPad], B::AlgMatElt[FldPad]) -> BoolElt, RngIntElt`

Return whether the axes of the hyperbolic isometries determined by the matrices intersect.
If so, then return the size of the path of the intersection.
Otherwise return the distance between the axes.

### Fixed sets on the Bruhat-Tits tree

`FixSet(tree::BTTree, A::AlgMatElt) -> BTTFix`

Give the fixed point set of the isometry on the tree.

The types of fixed points set are:
- `BTTFixSetIdentity` - The entire tree.
- `BTTFixSetBall` - A ball centred at a vertex.
- `BTTFixSetBallOnMidpoint` - A ball centered at the midpoint of an edge.
- `BTTFixSetHoroball`

`Center(fix::BTTFixSetBall) -> Tup`

`Center(fix::BTTFixSetBallOnMidpoint) -> Tup`

Return two vertices of maximal distance in the fixed ball.

`Radius(fix::BTTFixSetBall) -> RngIntElt`

`Radius(fix::BTTFixSetBallOnMidpoint) -> RngIntElt`

Return the radius of the fixed ball.

`BoundaryOnTree(fix::BTTFixSetHoroball) -> Tup`

Return a point on the tree contained in the boundary of the horoball.

`BoundaryAtInfinity(fix::BTTFixSetHoroball) -> Tup`

Return a point on the boundary of the tree fixed by the horoball.

`'*'(fix::BTTFixSetIdentity, mat::AlgMatElt[FldPad]) -> BTTFixSetIdentity`

`'*'(fix::BTTFixSetBand, mat::AlgMatElt[FldPad]) -> BTTFixSetBand`

`'*'(fix::BTTFixSetBall, mat::AlgMatElt[FldPad]) -> BTTFixSetBall`

`'*'(fix::BTTFixSetBallOnMidpoint, mat::AlgMatElt[FldPad]) -> BTTFixSetBallOnMidpoint`

`'*'(fix::BTTFixSetHoroball, mat::AlgMatElt[FldPad]) -> BTTFixSetHoroball`

Multiply the fixed point set by the matrix.
This gives the fixed point set of the conjugate of a corresponding isometry.

### Operations on boundary points

`CrossRatio(a::ModTupFldElt[FldPad], b::ModTupFldElt[FldPad], c::ModTupFldElt[FldPad], d::ModTupFldElt[FldPad]) -> FldPadElt`

Return the cross ratio between four points on the boundary of the tree.

`Midpoint(tree::BTTree, u::ModTupFldElt[FldPad], v::ModTupFldElt[FldPad], w::ModTupFldElt[FldPad]) -> BTTVertex`

Return the midpoint between 3 points on the boundary of the tree.

`ProjectionOntoMinTranslationSet(tree::BTTree, A::AlgMatElt[FldPad], v::ModTupFldElt[FldPad]) -> AlgMat[FldPad]`

Return the projection of a boundary point onto the minimum translation set of A.

### Other tree functions

`Random(tree::BTTree, radius::RngIntElt : boundary := false) -> BTTVert`

Return a random vertex of the tree with a maximum radius r.
If `boundary` is true, then return a random vertex with radius exactly r.