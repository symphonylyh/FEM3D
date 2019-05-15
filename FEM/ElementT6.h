// /**
//  * @file ElementT6.h
//  * Derived class from Element for the isoparametric T6 element.
//  *
//  * @author Haohang Huang
//  * @date Apr 22, 2018
//  */
//
// #ifndef ElementT6_h
// #define ElementT6_h
//
// #include "Element.h"
//
// /* Derived class for the isoparametric T6 element.
//  */
// class ElementT6 : public Element
// {
//     public:
//         /* See the documentation of base class Element. */
//         ElementT6();
//         ElementT6(const int & index, const std::vector<int> & nodeList, Node** const meshNode);
//         ~ElementT6();
//
//         Shape* const shape() const;
//
//     private:
//         /** A static structure that manages all the static members used in this class */
//         static staticMembers statics;
//
// };
//
// #endif /* ElementT6_h */
