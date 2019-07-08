#ifndef SET_H
#define SET_H

#include <Rcpp.h>
#include <vector>
#include <stack>
using namespace Rcpp;

// A container with smart memory allocation.
template <class T>
class Set {
        std::vector<T> items; // vector of items
        std::vector<int> active; // 1 or 0 to indicate whether each item is active or not
        std::stack<int> available; // ids of inactive items
        int n_active; // number of active items
        
        void validate(int id) {
            if ((id<0) || (id>=int(items.size())) || !active[id]) stop("(Set error) Given id does not exist");
        }
        

    public:
        static const int flag = -9; // special value to indicate conditions
        
        // Create new empty set
        Set() {
            n_active = 0;
        }

        int size() { return n_active; }
        
        // Return (pointer to) item with given id.
        T* item(int id) {
            validate(id);
            return &(items[id]); // & returns address of the item
        }
        
        // Activate a new element, and return the assigned id.
        int insert() {
            int id;
            if (!available.empty()) {
                id = available.top(); // get next available id from the stack
                available.pop(); // remove it from stack
                active[id] = 1; // indicate that this id is now active
            } else {
                // in this case, we must allocate more memory
                id = items.size();
                items.resize(id+1);
                active.push_back(1);
            }
            n_active += 1;
            return id;
        }
        
        // Remove element with given id.
        void remove(int id) {
            validate(id);
            available.push(id); // push the id onto the stack
            active[id] = 0; // indicate that this id is now inactive
            n_active -= 1; // decrement the count of active elements
        }
        
        // Return id of next active item after input id, or flag if there is no such item.
        int next(int id) {
            if (id<-1) stop("(Set error) Invalid input id for next()");
            while (1) {
                id += 1;
                if (id>=int(items.size())) return flag; // return flag if at the end of the list
                if (active[id]) return id; // otherwise, if at an active element, return the index
            }
        }
        
        // Return id of first active item.
        int first() { return next(-1); }
        
};


#endif













