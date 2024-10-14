#include <iostream>
#include <iomanip>
using namespace std;

bool isValidBase(char base)
{
    if (base == 'A' || base == 'T' || base == 'G' || base == 'C')
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool isValidStrand(string strand)
{
    if(strand == ""){
        return false;
    }
    for (int i = 0; i < (int)strand.length(); i++)
    {
        if (!(strand[i] == 'A' || strand[i] == 'T' || strand[i] == 'C' || strand[i] == 'G'))
        {
            return false;
        }
    }
    return true;
}


double strandSimilarity(string strand1, string strand2){
    int match = 0;
    if(strand1.length() != strand2.length()){
        return 0;
        
    }
    for(int i = 0; i < (int)strand1.length(); i++){
        if(strand1[i] == strand2[i]){
            match++;
            
        }
    }

    double similarity = (double)match/strand1.length();
    return similarity;
}


int bestStrandMatch(string input_strand,string target_strand, bool mutationAnalysis){
    if(input_strand.length()< target_strand.length()){

        cout << "Best similarity score: 0.0" << endl;  
        return -1;
    }
    double best_score = 0;
    int best_index = 0;

    for(int i = 0; i <= (int)input_strand.size()-(int)target_strand.size(); i ++){
         string sub_input = input_strand.substr(i, target_strand.size());
        
       
        double score = strandSimilarity(sub_input, target_strand);
        
      
        if (score > best_score) {
            best_score = score;
            best_index = i;
        }
    }

    
    cout << "Best similarity score: " << best_score << endl;
    if(mutationAnalysis){
    cout << "Best alignment index: " << best_index << endl;
    }
    return best_index;
    
}
void reverseComplement(string strand){
    string newStrand = "";
    //reverse strand
    for(int i = strand.length()-1; i >= 0; i --){
        newStrand += strand[i];
    }
    for(int i = 0;i < (int)strand.length(); i ++){
        if(newStrand[i]=='A'){
            newStrand[i]='T';
        }else if(newStrand[i]=='T'){
            newStrand[i]='A';
        }else if(newStrand[i]=='G'){
            newStrand[i]='C';
        }else if(newStrand[i]=='C'){
            newStrand[i]='G';
        }
    }
    cout << newStrand;

    


}

void transcribeDNAtoRNA(string strand){
    for(int i = 0; i < (int)strand.length(); i ++){
        if(strand[i]=='T'){
            strand[i] = 'U';
        }
    }
    cout << strand;
}

void identifyMutations(string input_strand, string target_strand) {
    int inputLength = input_strand.length();
    int targetLength = target_strand.length();
    bool swapped = false;
    // Swap if the input_strand is shorter than the target_strand
    if (inputLength < targetLength) {
        swap(input_strand, target_strand);
        swap(inputLength, targetLength);
        swapped = true;
    }

    int bestMatchIndex = bestStrandMatch(input_strand, target_strand,true);

    // Create a new target strand aligned with the best match
    string leftStars(bestMatchIndex, '*');
    string rightStars(inputLength - (bestMatchIndex + targetLength), '*');
    string alignedTargetStrand = leftStars + target_strand + rightStars;

    bool mutationsFound = false;
if(swapped==false){
    // Compare the aligned strands and detect mutations
    for (int i = 0; i < inputLength; i++) {
        if (alignedTargetStrand[i] == '*') {
            cout << "Deletion at position " << i+1 << ": " << input_strand[i] << " is deleted in target strand" << endl;
            mutationsFound = true;
        } else if (i >= targetLength + bestMatchIndex) {
            cout << "Insertion at position " << i+1 << ": " << alignedTargetStrand[i] << " is extra in the target strand" << endl;
            mutationsFound = true;
        } else if (input_strand[i] != alignedTargetStrand[i]) {
            cout << "Substitution at position " << i+1 << ": " << input_strand[i] << " -> " << alignedTargetStrand[i] << endl;
            mutationsFound = true;
        }
    }

    // If no mutations were found, print the message
    if (!mutationsFound) {
        cout << "No mutations found." << endl;
    }
}else if(swapped == true){
    for (int i = 0; i < inputLength; i++) {
        if (alignedTargetStrand[i] == '*') {
            cout << "Insertion at position " << i+1 << ": " << input_strand[i] << " is inserted in target strand" << endl;
            mutationsFound = true;
        } else if (i >= targetLength + bestMatchIndex) {
            cout << "Deletion at position " << i+1 << ": " << alignedTargetStrand[i] << " is deleted in the target strand" << endl;
            mutationsFound = true;
        } else if (input_strand[i] != alignedTargetStrand[i]) {
            cout << "Substitution at position " << i+1 << ": " << alignedTargetStrand[i] << " -> " << input_strand[i] << endl;
            mutationsFound = true;
        }
    }

    // If no mutations were found, print the message
    if (!mutationsFound) {
        cout << "No mutations found." << endl;
    }
}


}

void getCodingFrames(string strand) {
    int strandLength = strand.length();
    string startCodon = "ATG";
    string stopCodons[] = {"TAA", "TAG", "TGA"};
    bool foundFrame = false;

    for (int i = 0; i <= strandLength - 3; i++) {
       
        if (strand.substr(i, 3) == startCodon) {
            for (int j = i + 3; j <= strandLength - 3; j += 3) {
                
                string codon = strand.substr(j, 3);
                if (codon == stopCodons[0] || codon == stopCodons[1] || codon == stopCodons[2]) {
                    
                    if ((j - i) % 3 == 0) {
                        cout << strand.substr(i, j - i + 3) << endl;
                        foundFrame = true;
                        
                        i = j + 2; 
                        break; 
                    }
                }
            }
        }
    }

    if (!foundFrame) {
        cout << "The extracted reading frames are:" << endl;
        cout << "No reading frames found." << endl;
    }
}

void menu(){
    cout << "--- DNA Analysis Menu ---" << endl;
    cout << "1. Calculate the similarity between two sequences of the same length" << endl;
    cout << "2. Calculate the best similarity between two sequences of either equal or unequal length" << endl;
    cout << "3. Identify mutations" << endl;
    cout << "4. Transcribe DNA to RNA" << endl;
    cout << "5. Find the reverse complement of a DNA sequence" << endl;
    cout << "6. Extract coding frames" << endl;
    cout << "7. Exit" << endl;
    cout << "Please enter your choice (1 - 7):" << endl;

   
}

bool validSequence(string strand) {
    for (int i = 0; i < (int)strand.length(); i++) {
        if (!(strand[i] == 'A' || strand[i] == 'T' || strand[i] == 'G' || strand[i] == 'C')) {
            return false; 
        }
    }
    return true;  
}



int main() {
    int choice;
    string strand1, strand2;

    do {
        menu();
        cin >> choice;

        switch (choice) {
            case 1:
                do {
                    cout << "Enter the first DNA sequence: " << endl;
                    cin >> strand1;
                    if (!isValidStrand(strand1)) {
                        cout << "Invalid input. Please enter a valid sequence." << endl;
                    }
                } while (!isValidStrand(strand1));

                do {
                    cout << "Enter the second DNA sequence: " << endl;
                    cin >> strand2;
                    if (!isValidStrand(strand2)) {
                        cout << "Invalid input. Please enter a valid sequence." << endl;
                    }
                } while (!isValidStrand(strand2));

                if (strand1.length() != strand2.length()) {
                    cout << "Error: Input strands must be of the same length." << endl;
                } else {
                    cout << "Similarity score: " << strandSimilarity(strand1, strand2) << endl;
                }
                break;

            case 2:
                do {
                    cout << "Enter the first DNA sequence: " << endl;
                    cin >> strand1;
                    if (!isValidStrand(strand1)) {
                        cout << "Invalid input. Please enter a valid sequence." << endl;
                    }
                } while (!isValidStrand(strand1));

                do {
                    cout << "Enter the second DNA sequence: " << endl;
                    cin >> strand2;
                    if (!isValidStrand(strand2)) {
                        cout << "Invalid input. Please enter a valid sequence." << endl;
                    }
                } while (!isValidStrand(strand2));
                bestStrandMatch(strand1, strand2,false);
                
                
                break;

            case 3:
                do {
                    cout << "Enter the first DNA sequence: " << endl;
                    cin >> strand1;
                    if (!isValidStrand(strand1)) {
                        cout << "Invalid input. Please enter a valid sequence." << endl;
                    }
                } while (!isValidStrand(strand1));

                do {
                    cout << "Enter the second DNA sequence: " << endl;
                    cin >> strand2;
                    if (!isValidStrand(strand2)) {
                        cout << "Invalid input. Please enter a valid sequence." << endl;
                    }
                } while (!isValidStrand(strand2));
                identifyMutations(strand1, strand2);

                break;

            case 4:
                do {
                    cout << "Enter the DNA sequence to be transcribed:" << endl;
                    cin >> strand1;
                    if (!isValidStrand(strand1)) {
                        cout << "Invalid input. Please enter a valid sequence." << endl;
                    }
                } while (!isValidStrand(strand1));
                cout << "The transcribed DNA is: ";
                transcribeDNAtoRNA(strand1); 
                cout << endl; 
                break;

            case 5:
                do {
                    cout << "Enter the DNA sequence: " << endl;
                    cin >> strand1;
                    if (!isValidStrand(strand1)) {
                        cout << "Invalid input. Please enter a valid sequence." << endl;
                    }
                } while (!isValidStrand(strand1));
                cout << "The reverse complement is: ";
                reverseComplement(strand1);
                cout << endl;                    
                break;

            case 6:
                do {
                    cout << "Enter the DNA sequence: " << endl;
                    cin >> strand1;
                    if (!isValidStrand(strand1)) {
                        cout << "Invalid input. Please enter a valid sequence." << endl;
                    }
                } while (!isValidStrand(strand1));
                getCodingFrames(strand1);
                break;

            case 7:
                cout << "Exiting program." << endl;
                break;

            default:
                cout << "Invalid input. Please select a valid option." << endl;
                break;
        }
    } while (choice != 7);

    return 0;
}