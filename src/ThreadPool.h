#pragma once

#include <thread>
#include <mutex>
#include <queue>
#include <vector>
#include <functional>

class ThreadPool
{
public:
	void start();
	void QueueJob(std::function<void(int)> job, int id);
	void Stop();
	bool Busy();

	ThreadPool(int max_threads);
private:
	void ThreadLoop();

	bool Terminate = false;
	std::mutex m_Queue;
	std::condition_variable m_MutexCondition;
	std::vector<std::thread> m_threads;
	std::queue<std::pair<std::function<void(int)>, int>> jobs;

	int m_MaxThreads;

};

