#include "stepwidth_alignment_unit.hpp"



void ves::StepwidhtAlignmentUnit::setup(std::size_t interval, float min, float max, float target)
{
    tbb::spin_mutex::scoped_lock lock(mutex);
    alignment_every = interval;
    stepwidth = 0.5*(min+max);
    stepwidth_min = min;
    stepwidth_max = max;
    ratio_target = target;
    setup_flag = true;
}



float ves::StepwidhtAlignmentUnit::operator()() const
{
    return stepwidth;
}



void ves::StepwidhtAlignmentUnit::accepted()
{
    ++accepted_count;
    alignment_check();
}



void ves::StepwidhtAlignmentUnit::rejected()
{
    ++rejected_count;
    alignment_check();
}



float ves::StepwidhtAlignmentUnit::getTarget() const
{
    if( !setup_flag )
        throw std::logic_error("ves::StepwidhtAlignmentUnit::setup was not executed");
    else
        return ratio_target;
}



float ves::StepwidhtAlignmentUnit::getRatio() const
{   
    if( accepted_count == 0 && rejected_count == 0)
        return ratio_old;
    else
        return static_cast<float>(accepted_count) / (accepted_count+rejected_count);
}



std::size_t ves::StepwidhtAlignmentUnit::getAccepted() const
{
    return accepted_count;
}



std::size_t ves::StepwidhtAlignmentUnit::getRejected() const
{
    return rejected_count;
}



void ves::StepwidhtAlignmentUnit::alignment_check()
{
    if( !setup_flag )
        throw std::logic_error("ves::StepwidhtAlignmentUnit::setup was not executed");
    else if( accepted_count + rejected_count >= alignment_every ) 
    {
        // this will only be executed once per circle
        if(!mutex.try_lock())
            return;
        ratio_old = getRatio();
        do_aligment();
        accepted_count = rejected_count = 0;
        mutex.unlock();
    }
}


void ves::StepwidhtAlignmentUnit::do_aligment()
{
    float deviation = getRatio() - ratio_target;
    if( deviation < 0.0 )
    {
        if( stepwidth > stepwidth_min ) 
        {
            stepwidth *= (1.f+deviation);
        }
    }
    else
    {
        if( stepwidth >= stepwidth_max ) 
        {
            stepwidth = stepwidth_max;
        }
        else if( stepwidth < stepwidth_max )
        {
            stepwidth *= (1.f+deviation);
            if( stepwidth >= stepwidth_max ) 
            {
                stepwidth = stepwidth_max;
            }
        }
    }
}