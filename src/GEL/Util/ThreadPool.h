//
// Created by Cem Akarsubasi on 4/25/25.
//

#ifndef GEL_UTIL_THREADPOOL_H
#define GEL_UTIL_THREADPOOL_H

#include <thread>
#include <condition_variable>
#include <functional>
#include <queue>
#include <vector>

namespace GEL::Util
{

/// @brief a non-generic threadpool implementation
class ThreadPool
{
    std::mutex m_queue_mutex;
    std::condition_variable m_queue_condition;

    std::atomic_int m_number_working = 0;
    std::condition_variable m_number_working_condition;
    std::mutex m_waiting_mutex;

    std::queue<std::function<void()>> m_function_queue;
    std::vector<std::jthread> m_threads;

public:
    ThreadPool() = delete;

    /// A Threadpool should never be copied
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    // TODO: think about move construction

    /// @brief Construct a threadpool with the given thread count
    ///
    /// @throws std::invalid_argument if thread_count is 0
    explicit ThreadPool(const uint32_t thread_count) : m_threads{std::vector<std::jthread>(thread_count)}
    {
        if (thread_count == 0)
        {
            throw std::invalid_argument("thread_count must be greater than 0");
        }
        for (auto& thread : m_threads)
        {
            thread = std::jthread([this](const std::stop_token& stop_token)
            {
                while (!stop_token.stop_requested())
                {
                    std::function<void()> task;
                    {
                        std::unique_lock lock(m_queue_mutex);
                        m_queue_condition.wait(lock, [&]
                        {
                            return !m_function_queue.empty() || stop_token.stop_requested();
                        });

                        if (m_function_queue.empty())
                        {
                            continue;
                        }
                        task = std::move(m_function_queue.front());
                        m_function_queue.pop();
                    }
                    // TODO: exception handling to prevent a thread from going down
                    task();
                    m_number_working.fetch_sub(1, std::memory_order_release);
                    m_number_working_condition.notify_one();
                }
            });
        }
    }

    ~ThreadPool()
    {
        this->cancelAll();
        for (auto& thread : m_threads)
        {
            thread.join();
        }
    }

    /// @brief number of threads
    /// @return number of threads
    [[nodiscard]] size_t size() const
    {
        return m_threads.size();
    }

    /// @brief Adds a task to the queue
    /// @param task task to be executed
    /// @return the thread id
    void addTask(const std::function<void()>&& task)
    {
        m_number_working.fetch_add(1, std::memory_order_acquire);
        {
            std::lock_guard lock(m_queue_mutex);
            m_function_queue.push(task);
        }
        m_queue_condition.notify_one();
    }

    /// @brief Waits until all threads have finished
    void waitAll()
    {
        std::unique_lock lock(m_waiting_mutex);
        m_number_working_condition.wait(lock, [this]
        {
            return m_number_working.load(std::memory_order_relaxed) == 0;
        });
    }

    /// @brief Cancels every thread
    ///
    /// After this function, no other tasks should be added
    void cancelAll()
    {
        for (auto& t : m_threads)
        {
            t.request_stop();
        }
        m_queue_condition.notify_all();
    }
};

} // namespace GEL::Util

#endif // GEL_UTIL_THREADPOOL_H
