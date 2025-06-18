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
            workers.emplace_back([this] (int tid){
                //set_affinity(tid);
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
            },i);
        }
    }

    int set_affinity(int tid){
	    cpu_set_t mask;
	    CPU_ZERO(&mask);
	    CPU_SET(tid, &mask);
	    return sched_setaffinity(0, sizeof(mask), &mask);
    }

	size_t get_num_th() const{
			return this->num_threads;
	}

    template<class F, class... Args>
    void enfileira(F&& f, Args&&... args) {
       	auto task = std::make_shared<std::function<void()>>(std::bind(std::forward<F>(f), 
							    std::forward<Args>(args)...));
        {
		std::unique_lock<std::mutex> lock(queue_mutex);
		tasks.emplace([task](){ (*task)(); });
	}
       		condition.notify_one();
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

