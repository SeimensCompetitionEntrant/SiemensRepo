function complexConst = constManifoldMin(manifoldTruth,manifold) 
%Minimizes one manifold with a "truth" manifold up to a complex constant.
%Designed to take away the complex constant ambiguity inherent in our
%definition of a manifold

complexConst = (manifold'*manifold)^-1 * manifold'*manifoldTruth;