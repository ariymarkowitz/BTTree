Attach("tree.m");

Qp := pAdicField(2, 40);
Zp := Integers(Qp);
M := MatrixAlgebra(Qp, 2);

tree := BruhatTitsTree(Qp);

list := [M ! [Random(-10000, 10000) : i in [1..4]] : j in [1..10000]];
ellipticElements := [m : m in list | IsInvertible(m) and IsElliptic(m)];
balls := [m : m in ellipticElements | Type(FixSet(tree, m)) eq BTTFixSetBall or Type(FixSet(tree, m)) eq BTTFixSetBallOnMidpoint];
bands := [m : m in ellipticElements | Type(FixSet(tree, m)) eq BTTFixSetBand];

procedure timetest(l)
  a := 0;
  for m in l do
    a := FixSet(tree, m);
  end for;
end procedure;

t := Cputime();
timetest(ellipticElements);
print Cputime(t), Cputime(t)/#ellipticElements;

t := Cputime();
timetest(balls);
print Cputime(t), Cputime(t)/#balls;

t := Cputime();
timetest(bands);
print Cputime(t), Cputime(t)/#bands;
