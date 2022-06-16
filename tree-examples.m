Attach("tree.m");

Qp := pAdicField(2, 20);
M := MatrixAlgebra(Qp, 2);

tree := BruhatTitsTree(Qp);
print Origin(tree);

v := BTTVertexFromMatrix(tree, M ! [1, 2, 3, 22]);
print v;
print DistanceToOrigin(v);

A := M ! [2, 0, 1, 3];
print A*v;
print TranslationLength(A);

B := M ! [5, 3, 1, 3];
print TranslationLength(B);
print FixSet(tree, B);

print FixSet(tree, M ! [1, 0, 0, 1]);
print FixSet(tree, M ! [0, 1, 1, 0]);
print FixSet(tree, M ! [0, 1, -1, -2]);
print FixSet(tree, M ! [2, 1, 1/8, 5]);
print FixSet(tree, M ! [2, 1, 1/16, 5]);
print FixSet(tree, M ! [-3227, -4789, 4661, 7227]);