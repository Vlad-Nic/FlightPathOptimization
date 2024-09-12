#include <iostream>
#include <queue>
using namespace std;

// Function to find the kth largest element in an array using min heap
int findKthLargest(int arr[], int n, int k) {
    // Create a min heap of size k
    priority_queue<int, vector<int>, greater<int>> minHeap;

    // Insert the first k elements into the min heap
    for (int i = 0; i < k; ++i) {
        minHeap.push(arr[i]);
    }

    // For the remaining elements in the array
    for (int i = k; i < n; ++i) {
        // If the current element is greater than the root of the min heap
        // Replace the root with the current element and heapify
        if (arr[i] > minHeap.top()) {
            minHeap.pop();
            minHeap.push(arr[i]);
        }
    }

    // The root of the min heap will be the kth largest element
    return minHeap.top();
}

int main() {
    // Sample test cases
    int nums1[] = {3, 2, 1, 5, 6, 4};
    int k1 = 2;
    cout << "Output for [3, 2, 1, 5, 6, 4], k = 2: " << findKthLargest(nums1, 6, k1) << endl;

    int nums2[] = {3, 2, 3, 1, 2, 4, 5, 5, 6};
    int k2 = 4;
    cout << "Output for [3, 2, 3, 1, 2, 4, 5, 5, 6], k = 4: " << findKthLargest(nums2, 9, k2) << endl;

    int nums3[] = {10, 7, 11, 5, 27, 8, 9, 45};
    int k3 = 3;
    cout << "Output for [10, 7, 11, 5, 27, 8, 9, 45], k = 3: " << findKthLargest(nums3, 8, k3) << endl;

    int nums4[] = {100, 200, 300, 400, 500};
    int k4 = 2;
    cout << "Output for [100, 200, 300, 400, 500], k = 2: " << findKthLargest(nums4, 5, k4) << endl;

    int nums5[] = {1, 2, 3, 4, 5};
    int k5 = 5;
    cout << "Output for [1, 2, 3, 4, 5], k = 5: " << findKthLargest(nums5, 5, k5) << endl;

    return 0;
}
