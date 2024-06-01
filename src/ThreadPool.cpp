#include "ThreadPool.h"

void ThreadPool::start()
{
	for (int i = 0; i < m_MaxThreads; i++)
	{
		m_threads.emplace_back(std::thread(&ThreadPool::ThreadLoop, this));
	}
	Terminate = false;
}

void ThreadPool::QueueJob(std::function<void(int)> job, int id)
{
	{
		std::unique_lock<std::mutex> lock(m_Queue);
		jobs.push({ job,id });
	}
	m_MutexCondition.notify_one();
}

void ThreadPool::Stop()
{
	{
		std::unique_lock<std::mutex> lock(m_Queue);
		Terminate = true;
	}
	m_MutexCondition.notify_all();
	for (std::thread& active_thread : m_threads) {
		active_thread.join();
	}
	m_threads.clear();
}

bool ThreadPool::Busy()
{
	bool poolbusy;
	{
		std::unique_lock<std::mutex> lock(m_Queue);
		poolbusy = !jobs.empty();
	}
	return poolbusy;
}

ThreadPool::ThreadPool(int max_threads)
	:m_MaxThreads(max_threads)
{
}

void ThreadPool::ThreadLoop()
{
	while (true)
	{
		std::pair<std::function<void(int)>, int> job;
		{
			std::unique_lock<std::mutex> lock(m_Queue);
			m_MutexCondition.wait(lock, [this] {
				return !jobs.empty() || Terminate;
				});
			if (Terminate) {
				return;
			}
			job = jobs.front();
			jobs.pop();
		}
		job.first(job.second);
	}
}
