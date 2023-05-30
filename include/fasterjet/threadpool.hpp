/**
 * \file threadpool.hpp (description...)
 * 
 * \author A.S. Woodcock
 * 
 * \licence GPL3 (see COPYING.md)
*/

#pragma once
#include <thread>
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <functional>
#include <optional>
#include <semaphore>
#include <memory>
#include <deque>

using Task = std::function<void()>;

// a result shared between threads; automatically waits until avaliable before getting it
template <typename T>
class SharedResult
{
public:
   // set the return value and make it avaliable
   template<typename R>
   void set(R&& value)
   {
      // TODO: deal with already set
      result = std::forward<R>(value);
      ready.release();
   }
   // wait until return value is avaliable, then return it
   T get()
   {
      // TODO: deal with already get
      ready.acquire();
      return std::move(*result);
   }
private:
   std::binary_semaphore ready{ 0 };
   std::optional<T> result;
};

// specialization for void
template <>
class SharedResult<void>
{
public:
   void set() { ready.release(); }
   void get() { ready.acquire(); }
private:
   std::binary_semaphore ready{ 0 };
};

class ThreadPool
{
public:

   ThreadPool(const int64 workerCount)
   {
      workers.reserve(workerCount);
      for (int64 i=0; i<workerCount; ++i)
         workers.emplace_back(std::make_shared<std::jthread>(std::jthread(std::bind_front(&ThreadPool::threadKernel,this))));
   }

   template<typename F, typename...A>
   auto addTask(F&& function, A&&...args)
   {
      // the result will be stored in this shared pointer
      auto ret(std::make_shared<SharedResult<std::invoke_result_t<F, A...>>>());

      // wrap it up into a Task; movin in args by value
      auto task = [
         function = std::forward<F>(function),
         ...args = std::forward<A>(args),
         ret
      ]() mutable
      {
         if constexpr (std::is_void_v<std::invoke_result_t<F, A...>>)
         {
            function(std::forward<A>(args)...);
            ret->set();
         }
         else
         {
            ret->set(function(std::forward<A>(args)...));
         }
      };

      {
         std::lock_guard lk{ taskQueueMtx };
         tasks.push_back(std::move(task));
      }

      taskQueueCv.notify_one();
      return ret;
   }

   void waitForAllDone()
   {
      std::unique_lock lk{ taskQueueMtx };
      if (tasks.empty() && workingCount == 0) return;
      allDoneCv.wait(lk, [this] {return tasks.empty() && workingCount == 0; });
   }

   ~ThreadPool()
   {
      for (auto& w : workers)
         w->request_stop();
   }

private:

   void threadKernel(std::stop_token st)
   {
      while (true)
      {
         Task task;

         {
            std::unique_lock lk{ taskQueueMtx }; // locks the mutex
            if (tasks.empty()) // check if already task avaliable
               taskQueueCv.wait(lk, st, [this] { return !tasks.empty(); }); // unlock mutex and wait
            workingCount += 1;

            // check if we need to stop
            if (st.stop_requested())
               break;

            // get the task
            task = std::move(tasks.front());
            tasks.pop_front();
         }

         // do task
         task();
         workingCount -= 1;

         // // notify other threads when no tasks remaining
         // bool noTasks = false;
         // {
         //    std::unique_lock lk{ taskQueueMtx }; // locks the mutex
         //    noTasks = tasks.empty();
         // }

         {
            std::unique_lock lk{ taskQueueMtx }; // locks the mutex
            if (tasks.empty() && workingCount == 0)
               allDoneCv.notify_all();
         }
      }
   }

   // data
   std::atomic<int> workingCount = 0;
   std::mutex taskQueueMtx;
   // std::mutex allDoneMtx;
   std::condition_variable_any taskQueueCv;
   std::condition_variable allDoneCv;
   std::deque<Task> tasks;

   std::vector<std::shared_ptr<std::jthread>> workers;
};
