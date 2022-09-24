/**
* Utilities
*/

intrinsic EchelonBasis(V::ModTupRng) -> Tup
{ Give the basis of the vector space in reduced row echelon form }
  return Rows(EchelonForm(BasisMatrix(V)));
end intrinsic;

intrinsic EchelonForm(v::ModTupFldElt) -> ModTupFldElt
{ Return the reduced row echelon form of v }
  if v[1] eq 0 then
    return Parent(v) ! [0, 1];
  elif v[1] eq 1 then
    return v;
  else
    return v / v[1];
  end if;
end intrinsic;

intrinsic IsMultiple(v::ModTupFldElt, w::ModTupFldElt) -> BoolElt
{ Return true if w is a scalar multiple of v }
  return v[2]*w[1] eq v[1]*w[2];
end intrinsic;

intrinsic MinValuation(A::AlgMatElt[FldPad]) -> RngIntElt
{ Get the minimum valuation of a matrix over a p-adic field }
  return Min([Valuation(a) : a in ElementToSequence(A)]);
end intrinsic;

intrinsic VertexNormalForm(A::AlgMatElt[FldPad]) -> AlgMatElt[FldPad]
{ Convert a matrix to vertex normal form }
  K := BaseRing(A);
  p := UniformizingElement(K);
  n := MinValuation(A);
  B := ChangeRing(A * p^(-n), Integers(K));
  C := HermiteForm(B);
  return ChangeRing(C, K) * p^(-Valuation(C[1,1]));
end intrinsic;

intrinsic IsApproximatelyEqual(A::AlgMatElt[FldPad], B::AlgMatElt[FldPad]) -> BoolElt
{ Return true if A and B are equal, with different precisions }
  for i in [1..Rows(A)] do
    for j in [1..Columns(A)] do
      if not IsZero(A[i][j] - B[i][j]) then
        return false;
      end if;
    end for;
  end for;
  return true;
end intrinsic;

/**
* Bruhat-Tits tree
*/

declare type BTTree[BTTVert];
declare attributes BTTree: Field;
declare attributes BTTVert: Parent;
declare attributes BTTVert: Expansion;
declare attributes BTTVert: Precision;

intrinsic BruhatTitsTree(field::FldPad) -> BTTree
{ Create a Bruhat-Tits tree from a p-adic field }
  tree := New(BTTree);
  tree`Field := field;
  return tree;
end intrinsic;

intrinsic BruhatTitsTree(ring::RngPad) -> BTTree
{ Create a Bruhat-Tits tree from a p-adic ring }
  return BruhatTitsTree(FieldOfFractions(ring));
end intrinsic;

intrinsic Field(tree::BTTree) -> FldPad
{ Get the underlying field of the tree }
  return tree`Field;
end intrinsic;

intrinsic Prime(tree::BTTree) -> RngIntElt
{ Get a uniformizer of the tree }
  return UniformizingElement(Field(tree));
end intrinsic;

intrinsic Uniformizer(tree::BTTree) -> RngIntElt
{ Get a uniformizer of the tree }
  return UniformizingElement(Field(tree));
end intrinsic;

intrinsic Origin(tree::BTTree) -> BTTVert
{ Get the origin of the tree corresponding to the identity matrix }
  return BTTVertex(tree, 0, 0);
end intrinsic;

intrinsic 'eq'(x::BTTree, y::BTTree) -> BoolElt
{ Return true if the trees are equal }
  return x`Field eq y`Field;
end intrinsic;

intrinsic Print(T::BTTree)
{ Print T }
  printf "Bruhat-Tits Tree of %m", Field(T);
end intrinsic;

/**
* Bruhat-Tits tree vertex
*/

intrinsic BTTVertex(tree::BTTree, expansion::FldPadElt, precision::RngIntElt) -> BTTVert
{ Create a vertex as a p-adic approximation }
  return BTTVertex(tree, Matrix(Field(tree), [[1, expansion], [0, Prime(tree)^precision]]));
end intrinsic;

intrinsic BTTVertex(tree::BTTree, expansion::RngPadElt, precision::RngIntElt) -> BTTVert
{ Create a vertex as a p-adic approximation }
  return BTTVertex(tree, Matrix(Field(tree), [[1, expansion], [0, Prime(tree)^precision]]));
end intrinsic;

intrinsic BTTVertex(tree::BTTree, expansion::RngIntElt, precision::RngIntElt) -> BTTVert
{ Create a vertex as a p-adic approximation }
  return BTTVertex(tree, Matrix(Field(tree), [[1, expansion], [0, Prime(tree)^precision]]));
end intrinsic;

intrinsic BTTVertex(tree::BTTree, matrix::AlgMatElt[FldPad]) -> BTTVert
{ Convert a matrix to a Bruhat-Tits tree vertex }
  v := New(BTTVert);
  v`Parent := tree;
  mat := VertexNormalForm(matrix);
  error if mat[2][2] eq 0, "Matrix is singular";
  v`Expansion := mat[1][2];
  v`Precision := Valuation(mat[2][2]);
  return v;
end intrinsic;

intrinsic Matrix(v::BTTVert) -> AlgMatElt[FldPad]
{ Convert a vertex to a matrix in vertex normal form }
  return Matrix(Field(v), [[1, v`Expansion], [0, Prime(v)^(v`Precision)]]);
end intrinsic;

intrinsic Field(vertex::BTTVert) -> FldPad
{ Get the underlying field of a vertex }
  return Field(Parent(vertex));
end intrinsic;

intrinsic Prime(vertex::BTTVert) -> RngIntElt
{ Return the prime used for the valuation }
  return Prime(Parent(vertex));
end intrinsic;

intrinsic Uniformizer(vertex::BTTVert) -> RngIntElt
{ Return the prime used for the valuation }
  return Prime(Parent(vertex));
end intrinsic;

intrinsic Parent(vertex::BTTVert) -> BTTree
{ Get the Bruhat-Tits tree containing a vertex }
  return vertex`Parent;
end intrinsic;

intrinsic IsCoercible(T::BTTree, v::BTTVert) -> BoolElt, .
{ Return true if v is coercible to parent T with the conversion, otherwise return false }
  if T eq Parent(v) then
    return true, v;
  end if;
  F := Field(T);
  G := Field(v);
  b, u := IsCoercible(Expansion(v), G);
  if b then
    return true, BTTVertex(T, u, Precision(v));
  else
    return false;
  end if;
end intrinsic;

intrinsic '*'(v::BTTVert, mat::AlgMatElt) -> BTTVert
{ Multiply vertex by a matrix }
  return BTTVertex(Parent(v), Matrix(v)*mat);
end intrinsic;

intrinsic Expansion(v::BTTVert) -> FldPadElt
{ Return the expansion component of v }
  return v`Expansion;
end intrinsic;

intrinsic Precision(v::BTTVert) -> RngIntElt
{ Return the precision component of v }
  return v`Precision;
end intrinsic;

intrinsic 'eq'(x::BTTVert, y::BTTVert) -> BoolElt
{ Return true if the vertices are equal }
  return x`Parent eq y`Parent and IsZero(x`Expansion - y`Expansion) and x`Precision eq y`Precision;
end intrinsic;

intrinsic Print(v::BTTVert)
{ Print v }
  printf "Vertex %o mod p^%o", Expansion(v), Precision(v);
end intrinsic;

intrinsic Distance(v::BTTVert, w::BTTVert) -> RngIntElt
{ Return the distance between two vertices }
  return Precision(v) + Precision(w) - 2*Minimum([Valuation(Expansion(v)-Expansion(w)), Precision(v), Precision(w)]);
end intrinsic;

intrinsic Neighbor(v::BTTVert, x::ModTupFldElt[FldFin]) -> BTTVert
{ Return the neighbor w of v such that x is in the subspace v/w of k^2, where k is the residue field }
  error if IsZero(x), "Vector is zero";
  if x[1] eq 0 then
    return Neighbor(v, Infinity(x));
  else
    return Neighbor(v, x[2]/x[1]);
  end if;
end intrinsic;

intrinsic Neighbor(v::BTTVert, x::FldElt) -> BTTVert
{ Return the neighbor w of v such that [1, x] is in the subspace v/w of k^2, where k is the residue field }
  return BTTVertex(Parent(v), Expansion(v)+Field(v)!(x)*Prime(v)^Precision(v), Precision(v) + 1);
end intrinsic;

intrinsic Neighbor(v::BTTVert, x::RngElt) -> BTTVert
{ Return the neighbor w of v such that [1, x] is in the subspace v/w of k^2, where k is the residue field }
  return BTTVertex(Parent(v), Expansion(v)+Field(v)!(x)*Prime(v)^Precision(v), Precision(v) + 1);
end intrinsic;

intrinsic Neighbor(v::BTTVert, x::Infty) -> BTTVert
{ Return the neighbor w of v such that [0, 1] is in the subspace v/w of k^2, where k is the residue field }
  return BTTVertex(Parent(v), Expansion(v), Precision(v) - 1);
end intrinsic;

intrinsic TypeOfNeighbor(v::BTTVert, w::BTTVert) -> ModTupFldElt[FldFin]
{ Given neighbouring vertices v and w, reduced the vertex in v/w in Echelon form
  (ie. [1, u] or [0, 1], where u is in the residue field) }
  F := Field(v);
  p := Prime(v);
  k, m := ResidueClassField(Integers(Field(v)));
  V := VectorSpace(k, 2);
  if Precision(w) eq Precision(v) - 1 and Valuation(Expansion(v) - Expansion(w)) ge Precision(w) then
    return V![0, 1];
  else
    dif := Expansion(v) - Expansion(w);
    if Precision(w) eq Precision(v) + 1 and Valuation(dif) ge Precision(v) then
      return V![1, m(dif*Prime(v)^(-Precision(v)))];
    end if;
    error "w is not a neighbor of v";
  end if;
end intrinsic;

intrinsic Neighbors(v::BTTVert) -> SeqEnum[BTTVert]
{ Return the neighbours of v }
  k := ResidueClassField(Integers(Field(v)));
  return [Neighbor(v, x) : x in k] cat [Neighbor(v, Infinity())];
end intrinsic;

/**
* Properties of isometries
*/

intrinsic Mu(A::AlgMatElt) -> FldRatElt
{ Return Mu(A). Roughly speaking this is how far the characteristic polynomial of A is from being 'reduced'. }
  return Min(Valuation(Trace(A)), Valuation(Determinant(A))/2);
end intrinsic;

intrinsic TranslationLength(A::AlgMatElt) -> RngIntElt
{ Return the translation length of the action }
  require Nrows(A) eq 2 and Ncols(A) eq 2: "A is not a 2x2 matrix";
  return Integers() ! (Valuation(Determinant(A)) - 2*Mu(A));
end intrinsic;

intrinsic IsElliptic(mat::AlgMatElt) -> BoolElt
{ Return true if the isometry determined by mat is elliptic }
  return TranslationLength(mat) eq 0;
end intrinsic;

intrinsic IsHyperbolic(mat::AlgMatElt) -> BoolElt
{ Return true if the isometry determined by mat is hyperbolic }
  return TranslationLength(mat) ne 0;
end intrinsic;

/**
* Fixed sets on the Bruhat-Tits Tree
*/

declare type BTTFixSet;
declare attributes BTTFixSet: Tree;

declare type BTTFixSetIdentity: BTTFixSet;

declare type BTTFixSetBand: BTTFixSet;
declare attributes BTTFixSetBand: BoundaryPoints;
declare attributes BTTFixSetBand: Thickness;

declare type BTTFixSetBall: BTTFixSet;
declare attributes BTTFixSetBall: Center;
declare attributes BTTFixSetBall: Radius;

declare type BTTFixSetBallOnMidpoint: BTTFixSet;
declare attributes BTTFixSetBallOnMidpoint: Center;
declare attributes BTTFixSetBallOnMidpoint: Radius;

declare type BTTFixSetHoroball: BTTFixSet;
declare attributes BTTFixSetHoroball: BoundaryOnTree;
declare attributes BTTFixSetHoroball: BoundaryAtInfinity;

intrinsic Tree(fix::BTTFixSet) -> BTTree
{ Return the tree of the fixed point set }
  return fix`Tree;
end intrinsic;

intrinsic FixSetIdentity(tree::BTTree) -> BTTFixSetIdentity
{ Return a fixed point set of the Bruhat-Tits tree corresponding to the identity }
  fix := New(BTTFixSetIdentity);
  fix`Tree := tree;
  return fix;
end intrinsic;

intrinsic Print(x::BTTFixSetIdentity)
{ Print x }
  printf "Fixed set of identity";
end intrinsic;

intrinsic FixSetBand(tree::BTTree, nerveElt1::ModTupFldElt, nerveElt2::ModTupFldElt, thickness::RngIntElt) -> BTTFixSetBand
{ Return a fixed point set of the Bruhat-Tits tree corresponding to a band }
  fix := New(BTTFixSetBand);
  fix`Tree := tree;
  fix`BoundaryPoints := <nerveElt1, nerveElt2>;
  fix`Thickness := thickness;
  return fix;
end intrinsic;

intrinsic BoundaryPoints(fix::BTTFixSetBand) -> Tup
{ Return the points of the fixed band on the boundary at infinity, as vectors }
  return fix`BoundaryPoints;
end intrinsic;

intrinsic Thickness(fix::BTTFixSetBand) -> RngIntElt
{ Return the thickness of the fixed band }
  return fix`Thickness;
end intrinsic;

intrinsic Print(x::BTTFixSetBand)
{ Print x }
  printf "Fixed band with nerve between %o and thickness %o", BoundaryPoints(x), Thickness(x);
end intrinsic;

intrinsic FixSetBall(tree::BTTree, point::BTTVert, radius::RngIntElt) -> BTTFixSetBall
{ Return a fixed point set of the Bruhat-Tits tree corresponding to a ball with center at a vertex }
  fix := New(BTTFixSetBall);
  fix`Tree := tree;
  fix`Center := point;
  fix`Radius := radius;
  return fix;
end intrinsic;

intrinsic FixSetBallOnMidpoint(tree::BTTree, point1::BTTVert, point2::BTTVert, radius::FldRatElt) -> BTTFixSetBallOnMidpoint
{ Return a fixed point set of the Bruhat-Tits tree corresponding to a ball with center on the midpoint of an edge }
  fix := New(BTTFixSetBallOnMidpoint);
  fix`Tree := tree;
  fix`Center := <point1, point2>;
  assert Distance(point1, point2) eq 1;
  fix`Radius := radius;
  return fix;
end intrinsic;

intrinsic Center(fix::BTTFixSetBall) -> Tup
{ Return two vertices of maximal distance in the fixed ball }
  return fix`Center;
end intrinsic;

intrinsic Radius(fix::BTTFixSetBall) -> RngIntElt
{ Return the radius of the fixed ball }
  return fix`Radius;
end intrinsic;

intrinsic Print(x::BTTFixSetBall)
{ Print x }
  printf "Fixed ball with center %o and radius %o", Center(x), Radius(x);
end intrinsic;

intrinsic Center(fix::BTTFixSetBallOnMidpoint) -> Tup
{ Return two vertices of maximal distance in the fixed ball }
  return fix`Center;
end intrinsic;

intrinsic Radius(fix::BTTFixSetBallOnMidpoint) -> RngIntElt
{ Return the radius of the fixed ball }
  return fix`Radius;
end intrinsic;

intrinsic Print(x::BTTFixSetBallOnMidpoint)
{ Print x }
  printf "Fixed ball with center between %o and radius %o", Center(x), Radius(x);
end intrinsic;

intrinsic FixSetHoroball(tree::BTTree, vertex::BTTVert, boundary::ModTupFldElt) -> BTTFixSetHoroball
{ Return a fixed point set of the Bruhat-Tits tree corresponding to a horoball }
  fix := New(BTTFixSetHoroball);
  fix`Tree := tree;
  fix`BoundaryOnTree := vertex;
  fix`BoundaryAtInfinity := boundary;
  return fix;
end intrinsic;

intrinsic BoundaryOnTree(fix::BTTFixSetHoroball) -> Tup
{ Return a point on the tree contained in the boundary of the horoball. }
  return fix`BoundaryOnTree;
end intrinsic;

intrinsic BoundaryAtInfinity(fix::BTTFixSetHoroball) -> Tup
{ Return a point on the boundary of the tree fixed by the horoball. }
  return fix`BoundaryAtInfinity;
end intrinsic;

intrinsic Print(x::BTTFixSetHoroball)
{ Print x }
  printf "Fixed horoball bounded between %o and %o", BoundaryOnTree(x), BoundaryAtInfinity(x);
end intrinsic;

intrinsic '*'(fix::BTTFixSetIdentity, mat::AlgMatElt[FldPad]) -> BTTFixSetIdentity
{ Multiply the fixed point set by the matrix. This just returns the identity. }
  return fix;
end intrinsic;

intrinsic '*'(fix::BTTFixSetBand, mat::AlgMatElt[FldPad]) -> BTTFixSetBand
{ Multiply the fixed point set by the matrix.
This gives the fixed point set of the conjugate of a corresponding isometry. }
  p1 := BoundaryPoints(fix)[1];
  p2 := BoundaryPoints(fix)[2];
  return FixSetBand(Tree(fix), p1*mat, p2*mat, Thickness(fix));
end intrinsic;

intrinsic '*'(fix::BTTFixSetBall, mat::AlgMatElt[FldPad]) -> BTTFixSetBall
{ Multiply the fixed point set by the matrix.
This gives the fixed point set of the conjugate of a corresponding isometry. }
  return FixSetBall(Tree(fix), Center(fix)*mat, Radius(fix));
end intrinsic;

intrinsic '*'(fix::BTTFixSetBallOnMidpoint, mat::AlgMatElt[FldPad]) -> BTTFixSetBallOnMidpoint
{ Multiply the fixed point set by the matrix.
This gives the fixed point set of the conjugate of a corresponding isometry. }
  return FixSetBallOnMidpoint(Tree(fix), Center(fix)[1]*mat, Center(fix)[2]*mat, Radius(fix));
end intrinsic;

intrinsic '*'(fix::BTTFixSetHoroball, mat::AlgMatElt[FldPad]) -> BTTFixSetHoroball
{ Multiply the fixed point set by the matrix.
This gives the fixed point set of the conjugate of a corresponding isometry. }
  p1 := BoundaryOnTree(fix);
  p2 := BoundaryAtInfinity(fix);
  return FixSetHoroball(Tree(fix), p1*mat, p2*mat);
end intrinsic;

intrinsic FixSet(tree::BTTree, A::AlgMatElt) -> BTTFix
{ Give the fixed point set of the isometry on the tree }
  K := Field(tree);
  A := ChangeRing(A, K);

  assert IsElliptic(A);
  Z := Integers();

  if IsScalar(A) then
    return FixSetIdentity(tree);
  end if;

  KExt<a>, eigenvalues := SplittingField(CharacteristicPolynomial(A));
  eigenvectors := [EchelonBasis(Eigenspace(ChangeRing(A, KExt), x))[1] : x in eigenvalues];

  if KExt eq K then
    if #eigenvalues eq 2 then
      // There are 2 eigenvalues in K. The fixed set is a band.
      return FixSetBand(tree, eigenvectors[1], eigenvectors[2], Valuation(eigenvalues[1] - eigenvalues[2]) - Z!Mu(A));
    else
      // There is 1 eigenvalue in K. The fixed set is a horoball.
      // We need to find a point on the boundary of the fixed point set.
      // Project a vector v that is as far away from the eigenvalue as possible.
      v := Vector(K, (Valuation(eigenvectors[1][2]) gt 0) select [0, 1] else [1, 0]);
      return FixSetHoroball(tree, ProjectionOntoMinTranslationSet(tree, A, v), eigenvectors[1]);
    end if;
  else
    // (1, a) generates the ring of integers of KExt.
    // We may write a fixed boundary point x = alpha + beta*gen.
    // The lattice generated by {(1, x), (0, p'^n)} is a point in the original tree if and only if x = alpha mod (p'^n * O').
    e := RamificationIndex(KExt, K);
    fixedPointDecomp := Eltseq(eigenvectors[1][2]);
    n := Valuation(fixedPointDecomp[2])*e + Valuation(a);
    radius := (Valuation(eigenvalues[1] - eigenvalues[2]) - Valuation(eigenvectors[1][2] - eigenvectors[2][2]) + n)/e - Mu(A);

    if e eq 1 then
      return FixSetBall(tree, BTTVertex(tree, fixedPointDecomp[1], n), Z!radius);
    else
      v1 := BTTVertex(tree, fixedPointDecomp[1], n div 2);
      v2 := BTTVertex(tree, fixedPointDecomp[1], n div 2 + 1);
      return FixSetBallOnMidpoint(tree, v1, v2, radius);
    end if;
  end if;
end intrinsic;

/**
* Hyperbolic elements
*/

intrinsic TranslationAxisBoundary(A::AlgMatElt[FldPad]) -> Tup
{ Return the points on the boundary of the translation axis of the hyperbolic isometry determined by A }
  error if not IsHyperbolic(A), "A is not hyperbolic";
  K := BaseRing(A);
  roots := Roots(CharacteristicPolynomial(A));
  eig1 := EchelonBasis(Eigenspace(A, roots[1][1]))[1];
  eig2 := EchelonBasis(Eigenspace(A, roots[2][1]))[1];
  return <eig1, eig2>;
end intrinsic

intrinsic CrossRatio(a::ModTupFldElt[FldPad], b::ModTupFldElt[FldPad], c::ModTupFldElt[FldPad], d::ModTupFldElt[FldPad]) -> FldPadElt
{ Return the cross ratio between four points on the boundary of the tree }
  return ((a[2]*c[1] - a[1]*c[2])*(b[2]*d[1] - b[1]*d[2]))/((a[2]*b[1] - a[1]*b[2])*(c[2]*d[1] - c[1]*d[2]));
end intrinsic;

intrinsic Intersects(A::AlgMatElt[FldPad], B::AlgMatElt[FldPad]) -> BoolElt, RngIntElt
{ Return whether the axes of the hyperbolic isometries determined by the matrices intersect.
If so, returns the size of the path of the intersection.
Otherwise, returns the distance between the axes }
  points1 := TranslationAxisBoundary(A);
  points2 := TranslationAxisBoundary(B);
  cross1 := Valuation(CrossRatio(points1[1], points1[2], points2[1], points2[2]));
  cross2 := Valuation(CrossRatio(points1[2], points1[1], points2[1], points2[2]));

  m := Maximum(Abs(cross1), Abs(cross2));

  if (cross1 eq cross2) then
    if (cross1 eq 0) then
      return true, 0;
    else
      return false, m;
    end if;
  else
    return true, m;
  end if;
end intrinsic;

intrinsic CrossPath(tree::BTTree, A::AlgMatElt[FldPad], B::AlgMatElt[FldPad]) -> BoolElt, Tup
{ If the translation axes of A and B overlap, returns true and the endpoints of the overlapping path.
Otherwise, returns false and the endpoints of the minimal path between the axes. }
  points1 := TranslationAxisBoundary(A);
  points2 := TranslationAxisBoundary(B);
  a := points1[1];
  b := points1[2];
  c := points2[1];
  d := points2[2];
  aProj := Midpoint(tree, a, c, d);
  cProj := Midpoint(tree, c, a, b);
  if aProj eq cProj then
    bProj := Midpoint(tree, b, c, d);
    return true, <aProj, bProj>;
  else
    return false, <cProj, aProj>;
  end if;
end intrinsic;

/**
* Other vertex functions
*/

intrinsic IsometryBetweenAxes(tree::BTTree, v::ModTupFldElt[FldPad], w::ModTupFldElt[FldPad]) -> AlgMatElt[FldPad]
{ Return an isometry of translation length 1 inducing a path from v to w }
  v := EchelonForm(v);
  w := EchelonForm(w);
  K := Field(tree);
  D := Matrix(K, [[1, 0], [0, 2]]);
  error if v eq w, "v and w are scalar multiples of each other";
  M := Matrix([v, w]);
  return M^(-1)*D*M;
end intrinsic;

intrinsic IsometryBetweenVertices(v::BTTVert, w::BTTVert) -> SeqEnum[BTTVert]
{ Return an isometry of translation length 1 inducing a path from v to w }
  tree := Parent(v);
  K := Field(tree);
  A := IsometryBetweenAxes(tree, Vector(K, [1, Expansion(v)]), Vector(K, [1, Expansion(w)]));
  if Expansion(v) eq Expansion(w) and Precision(v) ge Precision(w) then
    A := A^-1;
  end if;
  return A;
end intrinsic;

intrinsic Path(v::BTTVert, w::BTTVert) -> SeqEnum[BTTVert]
{ Return the path between two vertices }
  error if Parent(v) ne Parent(w), "v and w have different parents";
  n := Minimum([Valuation(Expansion(v) - Expansion(w)), Precision(v), Precision(w)]);
  return [BTTVertex(Parent(v), Expansion(v), i) : i in [Precision(v)..n by -1]] cat [BTTVertex(Parent(w), Expansion(w), i) : i in [n+1..Precision(w)]];
end intrinsic

intrinsic IsInMinTranslationSet(A::AlgMatElt[FldPad], v::BTTVert) -> BoolElt
{ Return true if v is in the minimum translation set of A }
  return Distance(v, v*A) eq TranslationLength(A);
end intrinsic;

intrinsic IsInMinTranslationSetBoundary(A::AlgMatElt[FldPad], v::ModTupFldElt[FldPad]) -> BoolElt
{ Return true if v is in the minimum translation set of A }
  return IsMultiple(v*A, v);
end intrinsic;

intrinsic ProjectionOntoMinTranslationSet(tree::BTTree, A::AlgMatElt[FldPad], v::ModTupFldElt[FldPad]) -> AlgMat[FldPad]
{ Return the projection of v on the boundary onto A }
  error if IsInMinTranslationSetBoundary(A, v), "v is in the minimum translation set of A";
  p := Prime(tree);
  B := p^(Integers()!(-Mu(A))) * A;
  w := v * B;
  return BTTVertex(tree, Matrix([v, w]));
end intrinsic

intrinsic Midpoint(tree::BTTree, u::ModTupFldElt[FldPad], v::ModTupFldElt[FldPad], w::ModTupFldElt[FldPad]) -> BTTVertex
{ Return the midpoint between 3 points on the boundary of the tree }
  V := KernelMatrix(Matrix([u, v, w]))[1];
  return BTTVertex(tree, Matrix([u * V[1], v * V[2]]));
end intrinsic