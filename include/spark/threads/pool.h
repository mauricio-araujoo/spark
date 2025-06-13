#ifndef SPARK_THREADS_POOL_H
#define SPARK_THREADS_POOL_H

#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <future>
#include <memory>

namespace spark::threads {
class ThPool {
public:
    ThPool(size_t num_threads = std::thread::hardware_concurrency()) 
        : stop(false) {
		this->num_threads = num_threads;
        for(size_t i = 0; i < num_threads; ++i) {
            workers.emplace_back([this] {
                while(true) {
                    std::function<void()> task;
                    {
                        //TimeCounter t("TEMPO: ");
                        std::unique_lock<std::mutex> lock(queue_mutex);
                        condition.wait(lock, [this] {
                            return stop || !tasks.empty();
                        });
			            //delete t;
                        if(stop && tasks.empty())
                            return;

                        task = std::move(tasks.front());
                        tasks.pop();
                    }
                    task();
                }
            });
        }
    }
	size_t get_num_th() const{
			return this->num_threads;
	}

    template<class F, class... Args>
    auto enfileira(F&& f, Args&&... args) 
        -> std::future<typename std::result_of<F(Args...)>::type> {
        using return_type = typename std::result_of<F(Args...)>::type;
        
        auto task_ptr = std::make_shared<std::packaged_task<return_type()>>(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );
            
        std::future<return_type> res = task_ptr->get_future();
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            
            if(stop)
                throw std::runtime_error("enfileira em ThreadPool parada");
                
            tasks.emplace([task_ptr](){ (*task_ptr)(); });
        }
        condition.notify_one();
        return res;
    }

    ~ThPool() {
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            stop = true;
        }
        condition.notify_all();
        for(std::thread &worker: workers)
                worker.join();
    }

private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
    std::mutex queue_mutex;
    std::condition_variable condition;
	size_t num_threads;
    bool stop;
};
}

#endif

