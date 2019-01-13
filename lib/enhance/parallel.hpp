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

#include <tbb/tbb.h>
#include <memory>

#pragma GCC diagnostic ignored "-Weffc++" // hurts, but is necessary

namespace enhance
{
    // dummy for enqueue functions
    template<typename F>
    class lambda_task : public tbb::task 
    {
        F func_;
        tbb::task* execute() override { func_(); return NULL; }
    public:
        lambda_task( const F& __func) : func_(__func) { }
    };
    
    
    
    // enqueue root_task and use lambda function
    template<typename F>
    void enqueue_root( const F& __func ) 
    {
        tbb::task::enqueue( *new( tbb::task::allocate_root() ) enhance::lambda_task<F>(__func) );
    }
    
    
    
    // spawn root_task and use lambda function
    template<typename F>
    void spawn_root( const F& __func ) 
    {
        tbb::task::spawn( *new( tbb::task::allocate_root() ) enhance::lambda_task<F>(__func) );
    }
    
    
    
    // enqueue child_task and by using lambda function
    template<typename F>
    void enqueue_additional_child_of( const F& __func , tbb::empty_task& ROOT) 
    {
        tbb::task::enqueue( *new( tbb::task::allocate_additional_child_of(ROOT) ) enhance::lambda_task<F>(__func) );
    }
    
    
    
    // spawn child_task and by using lambda functionu
    template<typename F>
    void spawn_additional_child_of( const F& __func , tbb::empty_task& ROOT) 
    {
        tbb::task::spawn( *new( tbb::task::allocate_additional_child_of(ROOT) ) enhance::lambda_task<F>(__func) );
    }
    
    
    
    // spawns a root task with a given task F
    template<typename F = tbb::empty_task>
    struct root_task_base
    {
        F* operator->() { return task; }
        F& operator* () { return *task; }
        const F* operator->() const { return task; }
        const F& operator* () const { return *task; }
        
        template<typename...Args> inline void enqueue_child( Args...args ) { enhance::enqueue_additional_child_of(std::forward<Args>(args)..., *task); }
        template<typename...Args> inline void spawn_child  ( Args...args ) { enhance::spawn_additional_child_of  (std::forward<Args>(args)..., *task); }
        
    protected:
        root_task_base() : task ( new( tbb::task::allocate_root()) F ) {}
        virtual ~root_task_base(){ tbb::task::destroy(*task); } // destroy is critical to avoid memleak
        F* task;
    };
    
    
    // creates empty root task which will NOT be executed
    // can be waited for until all childs are finished
    // ###
    //   root_dummy ROOT;
    //   while(...)
    //   { 
    //       ROOT.enqueue_child([&]{FUNCTOR;}) 
    //   }
    //   wait_for_all() automatically when out of scope
    // ###
    struct root_dummy : public root_task_base<tbb::empty_task>
    {
        root_dummy() : root_task_base() { task->set_ref_count(1); }
    };
    
    
    
    struct scoped_root_dummy : public root_task_base<tbb::empty_task>
    {
        scoped_root_dummy() : root_task_base() { task->set_ref_count(1); }
        ~scoped_root_dummy(){ task->wait_for_all(); }
    };
    
    
    
    // spawn root task which will be executed
    template<typename F>
    struct root_task : public root_task_base<F>
    {
        root_task() : root_task_base<F>() { }
    };
    
    
    
    template<typename F>
    struct scoped_root_task : public root_task_base<F>
    {
        scoped_root_task() : root_task_base<F>() { }
        ~scoped_root_task(){ this->task->wait_for_all(); }
    };
}
