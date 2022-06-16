Attach("tree.m");

p := 2;
Qp := pAdicField(p, 20);
Zp := Integers(Qp);

for i in [1..12] do
  try
    P<x> := PolynomialAlgebra(Zp);
    Qp := SplittingField(P ! [Random(Zp) : i in [1..3]]);
    Zp := Integers(Qp);
    print(i);
  catch e;
    print "error at ", i;
    continue;
  end try;
end for;

M := MatrixAlgebra(Zp, 2);
T := BruhatTitsTree(Qp);

stop := false;
A := false;
while not stop do
  A := Random(M);
  if IsElliptic(A) then break; end if;
end while;

print(FixSet(T, A));
