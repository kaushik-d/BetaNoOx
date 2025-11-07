#ifndef beta_core_defines_h
#define beta_core_defines_h

class Node;
class BasicElement;
class Constraint;
class Load;
class NodeSet;
#include<vector>

typedef std::vector<Constraint*>	ConstraintList;
typedef std::vector<Load*>			LoadList;

#include "utility/array.hpp"
typedef Array<Node>					NodeGroup;
typedef Array<BasicElement>			ElementGroup;

typedef char ** charlist;

#define TEMPERATURE 0
#define MOISTURE    1

#define LINEAR 0
#define GEOMETRIC_NONLINEAR  1


#endif
