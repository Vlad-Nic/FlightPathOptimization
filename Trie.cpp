#include <iostream>
#include <string>
#include <vector>

// Trie node class
class TrieNode {
public:
    TrieNode() {
        for(int i = 0; i < ALPHABET_SIZE; i++) {
            children[i] = nullptr;
        }
        endOfWord = false;
    }

    static const int ALPHABET_SIZE = 26;
    TrieNode* children[ALPHABET_SIZE];
    bool endOfWord;
};

// Trie class
class Trie {
public:
    Trie(const std::vector<std::string>& words); // Constructor
    int count_words() const; // Method to count words
    void print() const; // Method to print words
    int longest_word() const; // Method to find the longest word

private:
    TrieNode* root;

    // Helper method to recursively count words
    int countWordsHelper(TrieNode* node) const;

    // Helper method to recursively print words
    void printHelper(TrieNode* node, std::string word) const;

    // Helper method to recursively find the longest word
    int longestWordHelper(TrieNode* node) const;

    // Helper method to insert a word into the Trie
    void insert(const std::string& word);
};

// Constructor implementation
Trie::Trie(const std::vector<std::string>& words): root(nullptr) {
    // Initialize the root node
    root = new TrieNode();

    // Insert each word from the vector into the Trie
    for (const std::string& word : words) {
        insert(word);
    }
}

// Method to insert a word into the Trie
void Trie::insert(const std::string& word) {
    TrieNode* cur = root;

    for (char ch : word) {
        int index = ch - 'a';
        if (!cur->children[index]) {
            cur->children[index] = new TrieNode();
        }
        cur = cur->children[index];
    }
    cur->endOfWord = true;
}

// Method to count words implementation
int Trie::count_words() const {
    // Start counting words from the root node
    return countWordsHelper(root);
}

// Helper method to recursively count words
int Trie::countWordsHelper(TrieNode* node) const {
    int count = 0;

    if (node == nullptr) return 0;

    // If current node is the end of a word, increment count
    if (node->endOfWord) count++;

    // Recursively count words in children nodes
    for (int i = 0; i < TrieNode::ALPHABET_SIZE; i++) {
        count += countWordsHelper(node->children[i]);
    }

    return count;
}

// Method to print words implementation
void Trie::print() const {
    // Start printing words from the root node
    printHelper(root, "");
}

// Helper method to recursively print words
void Trie::printHelper(TrieNode* node, std::string word) const {
    if (node == nullptr) return;

    // If current node is the end of a word, print it
    if (node->endOfWord) std::cout << word << std::endl;

    // Recursively print words in children nodes
    for (int i = 0; i < TrieNode::ALPHABET_SIZE; i++) {
        if (node->children[i]) {
            printHelper(node->children[i], word + char('a' + i));
        }
    }
}

// Method to find the longest word implementation
int Trie::longest_word() const {
    // Start finding the longest word from the root node
    return longestWordHelper(root);
}

// Helper method to recursively find the longest word
int Trie::longestWordHelper(TrieNode* node) const {
    int maxLength = 0;

    if (node == nullptr) return 0;

    // Check each child node for the longest word
    for (int i = 0; i < TrieNode::ALPHABET_SIZE; i++) {
        if (node->children[i]) {
            int length = 1 + longestWordHelper(node->children[i]);
            maxLength = std::max(maxLength, length);
        }
    }

    return maxLength;
}

int main() {
    std::vector<std::string> words = { "pickle", "squiggle", "bumblebee", "whippersnapper",
        "jellybean", "noodle", "kerplunk", "giggle",
        "muffin", "penguin"};
    Trie tr(words);
    std::cout << "There are " << tr.count_words() << " words in the Trie" << std::endl;
    std::cout << "Here are the words:\n";
    tr.print();
    std::cout << "The longest word has " << tr.longest_word() << " characters\n";
    return 0;
}
