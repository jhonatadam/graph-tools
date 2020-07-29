#ifndef GENERATOR_H
#define GENERATOR_H

#include <variant>
#include <functional>
#include <cassert>

using namespace std;

template<class Value, class State = Value>
class Generator
{
public:
    struct iterator {
        State state;
        function<Value(State const&)> value;
        function<void(State&)> next;
        function<bool(State const&, State const&)> comp;

        iterator& operator++()
        {
            next(state);
            return *this;
        }

        bool operator!=(iterator const& other) const
        {
            return !comp(state, other.state);
        }

        Value operator*() const
        {
            return value(state);
        }
    };

    Generator (State const& first, State const& last,
               function<Value (State const&)> value =
                       [](State const& s) -> Value {return s;},
               function<void(State &)> next =
                       [](State& s)->void{s++;},
               function<bool(State const&, State const&)> comp =
                       [](State const& s1, State const& s2) -> bool {return s1 == s2;})
{
    first_it.state = first;
    first_it.value = value;
    first_it.next = next;
    first_it.comp = comp;

    last_it.state = last;
}

iterator begin()
{
    return first_it;
}

iterator end()
{
    return last_it;
}

private:
iterator first_it;
iterator last_it;
};

#endif // GENERATOR_H
