/*  
*   Copyright 2019 Simon Raschke
*
*   Licensed under the Apache License, Version 2.0 (the "License");
*   you may not use this file except in compliance with the License.
*   You may obtain a copy of the License at
*
*       http://www.apache.org/licenses/LICENSE-2.0
*
*   Unless required by applicable law or agreed to in writing, software
*   distributed under the License is distributed on an "AS IS" BASIS,
*   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*   See the License for the specific language governing permissions and
*   limitations under the License.
*/

#pragma once

#include <vector>
#include <deque>
#include <list>
#include <algorithm>
#include <tbb/spin_rw_mutex.h>
#include <tbb/scalable_allocator.h>


namespace enhance
{
    // Class template for a thread safe container
    // tbb library required
    // locked by reader writer mutex
    // Notice: container will grow and shrink spin locked
    //         may result in overhead when accessed by many threads
    template
    <
        typename T, 
        template <typename, typename = tbb::scalable_allocator<T>> class Container
    >
    class ConcurrentContainer
    {
    public:
        typedef T member_t;
        typedef Container<T> container_t;
        typedef tbb::spin_rw_mutex mutex_t;
        typedef typename container_t::iterator iterator;
        typedef typename container_t::const_iterator const_iterator;

        ConcurrentContainer<T,Container>& operator=(ConcurrentContainer<T,Container>&&);

        // member management
        template<typename W>
        member_t& add(W&&);
        member_t& add();

        template<typename W>
        void add_if_new(W&&);
        template<typename W>
        void add_if_new(std::reference_wrapper<W>&&);

        template<typename W>
        bool contains(const W&);
        
        template<typename W>
        void remove(W&);
        void clear();

        // member access
        member_t& operator[](std::size_t);
        const member_t& operator[](std::size_t) const;
        member_t& back();
        const member_t& back() const;
        member_t& front();
        const member_t& front() const;
        container_t& data();
        const container_t& data() const;

        
        // information
        std::size_t size() const;
        bool empty() const;

        // iteration
        typename container_t::iterator begin();
        typename container_t::iterator end();
        typename container_t::const_iterator begin() const;
        typename container_t::const_iterator end() const;
        typename container_t::const_iterator cbegin() const;
        typename container_t::const_iterator cend() const;

        virtual ~ConcurrentContainer() = default;

    protected:
        container_t members {};

    private:
        // thread safety
        mutex_t mutex {};
    };



    template<typename T>
    using ConcurrentVector = ConcurrentContainer<T,std::vector>;
    template<typename T>
    using ConcurrentDeque = ConcurrentContainer<T,std::deque>;
    template<typename T>
    using ConcurrentList = ConcurrentContainer<T,std::list>;
    


    template<typename T,template <typename, typename> class Container>
    inline ConcurrentContainer<T,Container>& ConcurrentContainer<T,Container>::operator=(ConcurrentContainer<T,Container>&& other)
    {
        if(this!=&other)
        {
            members = other.members;
            other.members.clear();
        }
        return *this;
    }



    template<typename T,template <typename, typename> class Container>
    typename ConcurrentContainer<T,Container>::member_t& ConcurrentContainer<T,Container>::add()
    {
        mutex_t::scoped_lock lock(mutex, true);
        return members.emplace_back();
    }



    template<typename T,template <typename, typename> class Container>
    template<typename W>
    typename ConcurrentContainer<T,Container>::member_t& ConcurrentContainer<T,Container>::add(W&& other)
    {
        mutex_t::scoped_lock lock(mutex, true);
        return members.emplace_back(other);
    }



    template<typename T,template <typename, typename> class Container>
    template<typename W>
    void ConcurrentContainer<T,Container>::add_if_new(W&& other)
    {
        mutex_t::scoped_lock lock(mutex, false);
        if( !std::any_of(begin(),end(),[&](const T& member){ return member == other; }) )
        {   
            // vesDEBUG(__func__ << " other is new, so add")
            lock.upgrade_to_writer();
            members.emplace_back(other);
        }
        else
        {
            // vesDEBUG(__func__ << " other is not new, return");
        }
    }



    template<typename T,template <typename, typename> class Container>
    template<typename W>
    void ConcurrentContainer<T,Container>::add_if_new(std::reference_wrapper<W>&& other)
    {
        mutex_t::scoped_lock lock(mutex, false);
        if( !std::any_of(begin(),end(),[&](const T& member){ return member.get() == other.get(); }) )
        {   
            // vesDEBUG(__func__ << " other is new, so add")
            lock.upgrade_to_writer();
            members.emplace_back(other);
        }
        else
        {
            // vesDEBUG(__func__ << " other is not new, return");
        }
    }



    template<typename T,template <typename, typename> class Container>
    template<typename W>
    bool ConcurrentContainer<T,Container>::contains(const W& other)
    {
        mutex_t::scoped_lock lock(mutex, false);
        return std::any_of(begin(),end(),[&](const T& member){ return member == other; });
    }



    template<typename T,template <typename, typename> class Container>
    template<typename W>
    void ConcurrentContainer<T,Container>::remove(W& other)
    {
        mutex_t::scoped_lock lock(mutex, true);
        members.erase(std::remove(begin(), end(), other), end());
    }



    template<typename T,template <typename, typename> class Container>
    inline std::size_t ConcurrentContainer<T,Container>::size() const
    {
        return members.size();
    }



    template<typename T,template <typename, typename> class Container>
    inline bool ConcurrentContainer<T,Container>::empty() const
    {
        return members.empty();
    }



    template<typename T,template <typename, typename> class Container>
    inline void ConcurrentContainer<T,Container>::clear()
    {
        return members.clear();
    }



    template<typename T,template <typename, typename> class Container>
    inline typename ConcurrentContainer<T,Container>::member_t& ConcurrentContainer<T,Container>::operator[](std::size_t pos)
    {
        return members[pos];
    }



    template<typename T,template <typename, typename> class Container>
    inline const typename ConcurrentContainer<T,Container>::member_t& ConcurrentContainer<T,Container>::operator[](std::size_t pos) const
    {
        return members[pos];
    }



    template<typename T,template <typename, typename> class Container>
    inline typename ConcurrentContainer<T,Container>::member_t& ConcurrentContainer<T,Container>::back()
    {
        return members.back();
    }



    template<typename T,template <typename, typename> class Container>
    inline const typename ConcurrentContainer<T,Container>::member_t& ConcurrentContainer<T,Container>::back() const
    {
        return members.back();
    }



    template<typename T,template <typename, typename> class Container>
    inline typename ConcurrentContainer<T,Container>::member_t& ConcurrentContainer<T,Container>::front()
    {
        return members.front();
    }



    template<typename T,template <typename, typename> class Container>
    inline const typename ConcurrentContainer<T,Container>::member_t& ConcurrentContainer<T,Container>::front() const
    {
        return members.front();
    }



    template<typename T,template <typename, typename> class Container>
    inline typename ConcurrentContainer<T,Container>::container_t& ConcurrentContainer<T,Container>::data()
    {
        return members;
    }



    template<typename T,template <typename, typename> class Container>
    inline const typename ConcurrentContainer<T,Container>::container_t& ConcurrentContainer<T,Container>::data() const
    {
        return members;
    }



    template<typename T,template <typename, typename> class Container>
    inline typename ConcurrentContainer<T,Container>::container_t::iterator ConcurrentContainer<T,Container>::begin()
    {
        return std::begin(members);
    }



    template<typename T,template <typename, typename> class Container>
    inline typename ConcurrentContainer<T,Container>::container_t::iterator ConcurrentContainer<T,Container>::end()
    {
        return std::end(members);
    }



    template<typename T,template <typename, typename> class Container>
    inline typename ConcurrentContainer<T,Container>::container_t::const_iterator ConcurrentContainer<T,Container>::begin() const
    {
        return std::cbegin(members);
    }



    template<typename T,template <typename, typename> class Container>
    inline typename ConcurrentContainer<T,Container>::container_t::const_iterator ConcurrentContainer<T,Container>::end() const
    {
        return std::cend(members);
    }



    template<typename T,template <typename, typename> class Container>
    inline typename ConcurrentContainer<T,Container>::container_t::const_iterator ConcurrentContainer<T,Container>::cbegin() const
    {
        return std::cbegin(members);
    }



    template<typename T,template <typename, typename> class Container>
    inline typename ConcurrentContainer<T,Container>::container_t::const_iterator ConcurrentContainer<T,Container>::cend() const
    {
        return std::cend(members);
    }
}