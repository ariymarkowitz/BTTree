// TODO: In this file, matrices act via right multiplication.
// However, in Magma multiplication seems to be done on the left.

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

intrinsic MinValuation(A::AlgMatElt[FldPad]) -> RngIntEl
{ Get the minimum valuation of a matrix over a p-adic field }
  return Min([Valuation(a) : a in ElementToSequence(A)]);
end intrinsic;

intrinsic VertexNormalForm(A::AlgMatElt[FldPad]) -> AlgMatElt[FldPad]
{ Convert a matrix to vertex normal form }
  K := BaseRing(A);
  p := UniformizingElement(K);
  n := MinValuation(A);
  B := ChangeRing(A * p^(-n), Integers(K));
  C := Transpose(HermiteForm(Transpose(B)));
  return ChangeRing(C, K) * p^Valuation(C[1,1]);
end intrinsic;

declare type BTTree[BTTVert];
declare attributes BTTree: Field;
declare attributes BTTVert: Parent;
declare attributes BTTVert: Matrix;

/**
* Bruhat-Tits tree
*/

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

intrinsic Origin(tree::BTTree) -> BTTVert
{ Get the origin of the tree corresponding to the identity matrix }
  return BTTVertex(tree, 0, 0);
end intrinsic;

intrinsic 'eq'(x::BTTree, y::BTTree) -> BoolElt
{ Return true if the trees are equal }
  return x`Prime eq y`Prime;
end intrinsic;

intrinsic Print(T::BTTree)
{ Print T }
  printf "Bruhat-Tits Tree of %m", Field(T);
end intrinsic;

/**
* Bruhat-Tits tree vertex
*/

intrinsic BTTVertexFromMatrix(tree::BTTree, matrix::AlgMatElt[FldPad]) -> BTTVert
{ Convert a matrix to a Bruhat-Tits tree vertex }
  v := New(BTTVert);
  v`Parent := tree;
  mat := VertexNormalForm(matrix);
  error if mat[2][2] eq 0, "Matrix is singular";
  v`Matrix := VertexNormalForm(matrix);
  return v;
end intrinsic;

intrinsic BTTVertex(tree::BTTree, expansion::RngIntElt, precision::RngIntElt) -> BTTVert
{ Create a vertex as a p-adic approximation }
  return BTTVertexFromMatrix(tree, Matrix(Field(tree), [[1, 0], [expansion, Prime(tree)^precision]]));
end intrinsic;

intrinsic BTTVertex(tree::BTTree, expansion::FldPadElt, precision::RngIntElt) -> BTTVert
{ Create a vertex as a p-adic approximation }
  return BTTVertexFromMatrix(tree, Matrix(Field(tree), [[1, 0], [expansion, Prime(tree)^precision]]));
end intrinsic;

intrinsic Matrix(vertex::BTTVert) -> AlgMatElt[FldPad]
{ Convert a vertex to a matrix in vertex normal form }
  return vertex`Matrix;
end intrinsic;

intrinsic Field(vertex::BTTVert) -> FldPad
{ Get the underlying field of a vertex }
  return Field(Parent(vertex));
end intrinsic;

intrinsic Prime(vertex::BTTVert) -> RngIntElt
{ Return the prime used for the valuation }
  return Prime(Parent(vertex));
end intrinsic;

intrinsic Parent(vertex::BTTVert) -> BTTree
{ Get the Bruhat-Tits tree containing a vertex }
  return vertex`Parent;
end intrinsic;

intrinsic '*'(mat::AlgMatElt[FldPad], v::BTTVert) -> BTTVert
{ Multiply vertex by a matrix }
  return BTTVertexFromMatrix(Parent(v), mat*Matrix(v));
end intrinsic;

intrinsic Expansion(v::BTTVert) -> FldPadElt
{ Return the expansion component of v }
  return Matrix(v)[2, 1];
end intrinsic;

intrinsic Precision(v::BTTVert) -> RngIntElt
{ Return the precision component of v }
  return Valuation(Matrix(v)[2, 2]);
end intrinsic;

intrinsic 'eq'(x::BTTVert, y::BTTVert) -> BoolElt
{ Return true if the vertices are equal }
  return x`Parent eq y`Parent and x`Matrix eq y`Matrix;
end intrinsic;

intrinsic Print(v::BTTVert)
{ Print v }
  printf "Vertex %o mod p^%o", Expansion(v), Precision(v);
end intrinsic;

intrinsic DistanceToOrigin(v::BTTVert) -> RngaIntElt
{ Return the distance of a vertex from the origin }
  return Abs(Precision(v) - 2*Minimum(0, Valuation(Expansion(v))));
end intrinsic;

intrinsic Distance(v::BTTVert, w::BTTVert) -> RngIntElt
{ Return the distance between two vertices }
  return DistanceToOrigin(Matrix(w)^(-1) * v);
end intrinsic;

/**
* Properties of isometries
*/

intrinsic Mu(A::AlgMatElt) -> FldRatElt
{ Return Mu(A). Roughly speaking this is how far away the origin is from being on the boundary in a 'nice' way. }
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

intrinsic '*'(mat::AlgMatElt[FldPad], fix::BTTFixSetIdentity) -> BTTFixSetIdentity
{ Multiply the fixed point set by the matrix. This just returns the identity. }
  return fix;
end intrinsic;

intrinsic '*'(mat::AlgMatElt[FldPad], fix::BTTFixSetBand) -> BTTFixSetBand
{ Multiply the fixed point set by the matrix.
This gives the fixed point set of the conjugate of a corresponding isometry. }
  p1 := BoundaryPoints(fix)[1];
  p2 := BoundaryPoints(fix)[2];
  return FixSetBand(Tree(fix), p1*Transpose(mat), p2*Transpose(mat), Thickness(fix));
end intrinsic;

intrinsic '*'(mat::AlgMatElt[FldPad], fix::BTTFixSetBall) -> BTTFixSetBall
{ Multiply the fixed point set by the matrix.
This gives the fixed point set of the conjugate of a corresponding isometry. }
  return FixSetBall(Tree(fix), mat*Center(fix), Radius(fix));
end intrinsic;

intrinsic '*'(mat::AlgMatElt[FldPad], fix::BTTFixSetBallOnMidpoint) -> BTTFixSetBallOnMidpoint
{ Multiply the fixed point set by the matrix.
This gives the fixed point set of the conjugate of a corresponding isometry. }
  return FixSetBallOnMidpoint(Tree(fix), mat*Center(fix)[1], mat*Center(fix)[2], Radius(fix));
end intrinsic;

intrinsic '*'(mat::AlgMatElt[FldPad], fix::BTTFixSetHoroball) -> BTTFixSetHoroball
{ Multiply the fixed point set by the matrix.
This gives the fixed point set of the conjugate of a corresponding isometry. }
  p1 := BoundaryOnTree(fix);
  p2 := BoundaryAtInfinity(fix);
  return FixSetHoroball(Tree(fix), mat*p1, p2*Transpose(mat));
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
      p := UniformizingElement(K);
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
  return <Vector(K, [1, roots[1][1]]), Vector(K, [1, roots[2][1]])>;
end intrinsic

intrinsic CrossRatio(a::ModTupFldElt[FldPad], b::ModTupFldElt[FldPad], c::ModTupFldElt[FldPad], d::ModTupFldElt[FldPad]) -> FldPadElt
{ Return the cross ratio between four points on the boundary of the tree }
  return ((a[2]*c[1] - a[1]*c[2])*(b[2]*d[1] - b[1]*d[2]))/((a[2]*b[1] - a[1]*b[2])*(c[2]*d[1] - c[1]*d[2]));
end intrinsic;

intrinsic Intersects(A::AlgMatElt[FldPad], B::AlgMatElt[FldPad]) -> BoolElt, RngIntElt
{ Return whether the axes of the isometries determined by the matrices intersect. If so, returns the size of the path of the intersection.
Otherwise, returns the distance between the axes }
  points1 := TranslationAxisBoundary(A);
  points2 := TranslationAxisBoundary(B);
  cross1 := Valuation(CrossRatio(points1[1], points1[2], points2[1], points2[2]));
  cross2 := Valuation(CrossRatio(points1[2], points1[1], points2[1], points2[2]));

  m := Maximum(cross1, cross2);

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
Otherwise, returns false and the endpoints of the minimal path between the axes.}
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
    return false, <aProj, cProj>;
  end if;
end intrinsic;

/**
* Other vertex functions
*/

intrinsic IsometryBetweenAxes(tree::BTTree, v::ModTupFldElt[FldPad], w::ModTupFldElt[FldPad]) -> AlgMatElt[FldPad]
{ Return an isometry of translation length 2 inducing a path from v to w }
  v := EchelonForm(v);
  w := EchelonForm(w);
  K := Field(tree);
  D := Matrix(K, [[1, 0], [0, 2]]);
  error if v eq w, "v and w are scalar multiples of each other";
  M := Transpose(Matrix(K, [v, w]));
  return M*D*M^(-1);
end intrinsic;

intrinsic IsometryBetweenVertices(v::BTTVert, w::BTTVert) -> SeqEnum[BTTVert]
{ Return an isometry of translation length 2 inducing a path from v to w }
  error if Parent(v) ne Parent(w), "v and w do not have the same parent";
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
  A := IsometryBetweenVertices(v, w);
  return [A^i * v : i in [1 .. Distance(v, w)]];
end intrinsic

intrinsic IsInMinTranslationSet(A::AlgMatElt[FldPad], v::BTTVert) -> BoolElt
{ Return true if v is in the minimum translation set of A }
  return Distance(v, A*v) eq TranslationLength(A);
end intrinsic;

intrinsic IsInMinTranslationSetBoundary(A::AlgMatElt[FldPad], v::ModTupFldElt[FldPad]) -> BoolElt
{ Return true if v is in the minimum translation set of A }
  return IsMultiple(v*Transpose(A), v);
end intrinsic;

intrinsic ProjectionOntoMinTranslationSet(tree::BTTree, A::AlgMatElt[FldPad], v::ModTupFldElt[FldPad]) -> AlgMat[FldPad]
{ Return the projection of v on the boundary onto A }
  error if IsInMinTranslationSetBoundary(A, v), "v is in the minimum translation set of A";
  K := Field(tree);
  p := UniformizingElement(K);
  B := A * p^(Integers()!(-Mu(A)));
  w := v * Transpose(B);
  return BTTVertexFromMatrix(tree, Transpose(Matrix([v, w])));
end intrinsic

intrinsic Midpoint(tree::BTTree, u::ModTupFldElt[FldPad], v::ModTupFldElt[FldPad], w::ModTupFldElt[FldPad]) -> BTTVertex
{ Return the midpoint between 3 points on the boundary of the tree }
  V := KernelMatrix(Matrix([u, v, w]))[1];
  return BTTVertexFromMatrix(tree, Transpose(Matrix(V[1] * u, V[2] * v)));
end intrinsic