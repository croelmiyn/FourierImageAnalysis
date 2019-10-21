//////////////////////////////////////////////////////////////////////////////////
// Concurrency utilities														//
// Bundled from code written by Piotr Wendykier      							//
// CC BY SA	by Remy Colin and the Max Planck Society      						//
//////////////////////////////////////////////////////////////////////////////////
// Date   : 2019-10-21															//
// Author : Rémy Colin															//
// Email  : remycolin@laposte.net / remy.colin@synmikro.mpi-marburg.mpg.de		//
// Adress : Max Planck Institute for Terrestrial Microbiology					//
// 			Karl-von-Frisch-strasse 10											//
//          35043 Marburg, Germany												//
//////////////////////////////////////////////////////////////////////////////////

package mpi.rc.IJ.IJutilities;

import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.TimeUnit;

/**
 * Concurrency utilities.
 * <p>
 * @author Piotr Wendykier (p.wendykier@icm.edu.pl)
 */
public class ConcurrencyUtils
{

    /**
     * Thread pool.
     */
    private static final ExecutorService DEFAULT_THREAD_POOL = Executors.newCachedThreadPool(new ConcurrencyUtils.CustomThreadFactory(new ConcurrencyUtils.CustomExceptionHandler()));

    private static ExecutorService threadPool = DEFAULT_THREAD_POOL;

    private static int nthreads = getNumberOfProcessors();

    private static long concurrentThreshold = 100000;

    private ConcurrencyUtils()
    {

    }

    private static class CustomExceptionHandler implements Thread.UncaughtExceptionHandler
    {

        @Override
        public void uncaughtException(Thread t, Throwable e)
        {
            e.printStackTrace();
        }

    }

    private static class CustomThreadFactory implements ThreadFactory
    {

        private static final ThreadFactory DEFAULT_FACTORY = Executors.defaultThreadFactory();

        private final Thread.UncaughtExceptionHandler handler;

        CustomThreadFactory(Thread.UncaughtExceptionHandler handler)
        {
            this.handler = handler;
        }

        @Override
        public Thread newThread(Runnable r)
        {
            Thread t = DEFAULT_FACTORY.newThread(r);
            t.setUncaughtExceptionHandler(handler);
            return t;
        }
    };

    /**
     * Returns the minimum length of array for which multiple threads are used.
     * <p>
     * @return the minimum length of array for which multiple threads are used
     */
    public static long getConcurrentThreshold()
    {
        return ConcurrencyUtils.concurrentThreshold;
    }

    /**
     * Sets the minimum length of an array for which multiple threads are used.
     * <p>
     * @param concurrentThreshold minimum length of an array for which multiple threads are used
     */
    public static void setConcurrentThreshold(long concurrentThreshold)
    {
        ConcurrencyUtils.concurrentThreshold = max(1, concurrentThreshold);
    }

    /**
     * Returns the number of available processors.
     * <p>
     * @return number of available processors
     */
    public static int getNumberOfProcessors()
    {
        return Runtime.getRuntime().availableProcessors();
    }

    /**
     * Returns the current number of threads.
     * <p>
     * @return the current number of threads.
     */
    public static int getNumberOfThreads()
    {
        return ConcurrencyUtils.nthreads;
    }

    /**
     * Sets the number of threads.
     * <p>
     * @param n new value of threads
     */
    public static void setNumberOfThreads(int n)
    {
        ConcurrencyUtils.nthreads = n;
    }

    /**
     * Submits a value-returning task for execution and returns a Future
     * representing the pending results of the task.
     * <p>
     * @param <T>  type
     * @param task task for execution
     * <p>
     * @return handle to the task submitted for execution
     */
    public static <T> Future<T> submit(Callable<T> task)
    {
        if (ConcurrencyUtils.threadPool.isShutdown() || ConcurrencyUtils.threadPool.isTerminated()) {
            ConcurrencyUtils.threadPool = DEFAULT_THREAD_POOL;
        }
        return ConcurrencyUtils.threadPool.submit(task);
    }

    /**
     * Submits a Runnable task for execution and returns a Future representing that task.
     * <p>
     * @param task task for execution
     * <p>
     * @return handle to the task submitted for execution
     */
    public static Future<?> submit(Runnable task)
    {
        if (ConcurrencyUtils.threadPool.isShutdown() || ConcurrencyUtils.threadPool.isTerminated()) {
            ConcurrencyUtils.threadPool = DEFAULT_THREAD_POOL;
        }
        return ConcurrencyUtils.threadPool.submit(task);
    }

    /**
     * Waits for all threads to complete computation.
     * <p>
     * @param futures list of handles to the tasks
     * <p>
     * @throws ExecutionException   if the computation threw an exception
     * @throws InterruptedException if the current thread was interrupted while waiting
     */
    public static void waitForCompletion(Future<?>[] futures) throws InterruptedException, ExecutionException
    {
        int size = futures.length;
        for (int j = 0; j < size; j++) {
            futures[j].get();
        }
    }

    /**
     * Sets the pool of threads.
     * <p>
     * @param threadPool pool of threads
     */
    public static void setThreadPool(ExecutorService threadPool)
    {
        ConcurrencyUtils.threadPool = threadPool;
    }

    /**
     * Returns the pool of threads.
     * <p>
     * @return pool of threads
     */
    public static ExecutorService getThreadPool()
    {
        return ConcurrencyUtils.threadPool;
    }

    /**
     * Shutdowns all submitted tasks.
     */
    public static void shutdownThreadPoolAndAwaitTermination()
    {
        ConcurrencyUtils.threadPool.shutdown(); // Disable new tasks from being submitted
        try {
            // Wait a while for existing tasks to terminate
            if (!ConcurrencyUtils.threadPool.awaitTermination(60, TimeUnit.SECONDS)) {
                ConcurrencyUtils.threadPool.shutdownNow(); // Cancel currently executing tasks
                // Wait a while for tasks to respond to being cancelled
                if (!ConcurrencyUtils.threadPool.awaitTermination(60, TimeUnit.SECONDS))
                    System.err.println("Pool did not terminate");
            }
        } catch (InterruptedException ie) {
            // (Re-)Cancel if current thread also interrupted
            ConcurrencyUtils.threadPool.shutdownNow();
            // Preserve interrupt status
            Thread.currentThread().interrupt();
        }
    }

    private static long max(long a, long b){
        return (a>b)?a:b;
    }
}

