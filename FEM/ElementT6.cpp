// /**
//  * @file ElementT6.cpp
//  * Implementation of ElementT6 class.
//  *
//  * @author Haohang Huang
//  * @date Apr 22, 2018
//  */
//
// #include <stdio.h>
//
// #include "ElementT6.h"
//
// ElementT6::staticMembers ElementT6::statics(8,9,4,3,3);
//
// ElementT6::ElementT6()
// {
// }
//
// ElementT6::ElementT6(const int & index, const std::vector<int> & nodeList, Node** const meshNode)
//     : Element(index, nodeList, meshNode) // call the constructor of base class in the initializer list!
// {
// }
//
// ElementT6::~ElementT6()
// {
// }
//
// Shape* const ElementT6::shape() const
// {
//     return statics.shape;
// }
