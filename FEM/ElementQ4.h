/**
 * @file ElementQ4.h
 * Derived class from Element for the isoparametric Q4 element.
 *
 * @author Haohang Huang
 * @date Feburary 10, 2018
 * @note Efficiency optimized by polymorph shape on March 26, 2018.
 * @note Efficiency optimized by storing local stiffness matrix and return-by-ref
 * on March 27, 2018
 * @note Efficiency optimized by the generalization of all element-wise operations
 * into base class Element on Apr 22, 2018.
 */

#ifndef ElementQ4_h
#define ElementQ4_h

#include "Element.h"

/* Derived class for the isoparametric Q4 element.
 */
class ElementQ4 : public Element
{
    public:
        /* See the documentation of base class Shape. */
        ElementQ4();
        ElementQ4(const int & index, const std::vector<int> & nodeList, Node** const meshNode, Material* const material);
        ~ElementQ4();

        Shape* shape() const;

    private:
        /** A static structure that manages all the static members used in this class */
        static staticMembers statics;

};

#endif /* ElementQ4_h */
