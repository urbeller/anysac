# anysac
AnySac: a collection of Ransac like algorithms.

For this first version, I have only implemented the classical [Ransac](https://en.wikipedia.org/wiki/Random_sample_consensus). To use it, user must supply a driver object that must implement these methods (comments show an example for line fitting):

```c++
 int model_size() const; // For line fitting returns 2.
 int data_size() const; // Returns the data size.
 bool fit_model(const vector<int> &indices, vector<ModelType> &models); // This takes indices and returns the possible lines.
 int count_inliers(const ModelType &model); // The number of data points that agree with this model.
 bool check_subset(const vector<int> &indices); // Is the subset good ? ex: make sure that 2 points ar different.
 bool has_generator(); // if true, Ransac will call the supplied get_subset().
 bool get_subset(vector<int> &indices); // Gives a random combination of indices. For lines, it will return 2 random indices.
```
Please have a look at the [line fitting](src/line_driver.cxx) example.
