#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/mman.h>
#include <unistd.h>
#include <time.h>

typedef unsigned short WORD;
typedef unsigned int DWORD;
typedef unsigned int LONG;
typedef unsigned char byte;

#pragma pack(push, 1)
typedef struct tagBITMAPFILEHEADER
{
    WORD bfType;      //specifies the file type
    DWORD bfSize;     //specifies the size in bytes of the bitmap file
    WORD bfReserved1; //reserved; must be 0
    WORD bfReserved2; //reserved; must be 0
    DWORD bfOffBits;  //species the offset in bytes from the bitmapfileheader to the bitmap bits
} tagBITMAPFILEHEADER;
#pragma pack(pop)

#pragma pack(push, 1)
typedef struct tagBITMAPINFOHEADER
{
    DWORD biSize;         //specifies the number of bytes required by the struct
    LONG biWidth;         //specifies width in pixels
    LONG biHeight;        //species height in pixels
    WORD biPlanes;        //specifies the number of color planes, must be 1
    WORD biBitCount;      //specifies the number of bit per pixel
    DWORD biCompression;  //spcifies the type of compression
    DWORD biSizeImage;    //size of image in bytes
    LONG biXPelsPerMeter; //number of pixels per meter in x axis
    LONG biYPelsPerMeter; //number of pixels per meter in y axis
    DWORD biClrUsed;      //number of colors used by th ebitmap
    DWORD biClrImportant; //number of colors that are important
} tagBITMAPINFOHEADER;
#pragma pack(pop)

typedef struct Node
{
    int pixel;
    float freq;
    struct Node *left, *right;
    char huffman_code[100];
} Node;

typedef struct HuffTree
{
    int pixel, arr_loc;
    float freq;
} HuffTree;

FILE *load_header(char *filename, tagBITMAPINFOHEADER *bitmapInfoHeader, tagBITMAPFILEHEADER *bitmapFileHeader) {
    FILE *file; //our file pointer
    unsigned char *bmp_image; //store image data

    //open filename in read binary mod
    file = fopen(filename, "rb");
    if (file == NULL)
    {
        printf("File do not exist!\n");
        exit(1);
    }

    //read the bitmap file header
    fread(bitmapFileHeader, sizeof(tagBITMAPFILEHEADER), 1, file);

    //verify that this is a bmp file by check bitmap id
    if (bitmapFileHeader->bfType != 0x4D42)
    {
        fclose(file);
        printf("File is not a .bmp type!\n");
        exit(1);
    }

    //read the bitmap info header
    fread(bitmapInfoHeader, sizeof(tagBITMAPINFOHEADER), 1, file);

    return file;
}

byte *load_bmp_file(char *filename, tagBITMAPINFOHEADER *bitmapInfoHeader, tagBITMAPFILEHEADER *bitmapFileHeader)
{
    FILE *file = load_header(filename, bitmapInfoHeader, bitmapFileHeader); //our file pointer
    unsigned char *bmp_image; //store image data

    //move file point to the begging of bitmap data
    fseek(file, bitmapFileHeader->bfOffBits, SEEK_SET);

    //allocate enough memory for he bitmap image data
    bmp_image = (unsigned char *)mmap(NULL, bitmapInfoHeader->biSizeImage, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

    //read in the bitmap image data
    fread(bmp_image, bitmapInfoHeader->biSizeImage, 1, file);

    //make sure bitmap image data was read
    if (bmp_image == NULL)
    {
        fclose(file);
        return NULL;
    }

    //close file and return bitmap iamge data
    fclose(file);

    return bmp_image;
}

void verify_input(int argc, char *argv[], int *format)
{
    if (argc != 4) {
        printf("Please input 4 arguments! Please input as format like below (extension: .bmp to compress OR .hbmp to decompress):\nmycompress [flags] [filename] [encryption password]\n");
        exit(1);
    }

    if (strcmp(argv[1], "g") != 0 && strcmp(argv[1], "c") != 0) {
        printf("Please input the correct flags for the colors! 'g' for grey or 'c' for color.\n");
        exit(1);
    }

    char *pic = (char *)malloc(strlen(argv[2]) + 1);
    strcpy(pic, argv[2]);

    char *ptr_pic_bmp = strstr(pic, ".bmp"), *ptr_pic_hbmp = strstr(pic, ".hbmp");

    if (ptr_pic_bmp == NULL && ptr_pic_hbmp == NULL) {
        printf("Invalid input! Please input as format like below (extension: .bmp to compress OR .hbmp to decompress):\nmycompress [flags] [filename] [encryption password]\n");
        exit(1);
    }
    if (ptr_pic_bmp != NULL)
        *format = 0; //user wanted to compress

    if (ptr_pic_hbmp != NULL)
        *format = 1; //user wanted to decompress

    if (strlen(argv[3]) != 7) {
        printf("Password must be exactly 7 characters!\n");
        exit(1);
    }

    free(pic);
}

void encrypt_and_decrypt(char *pass, char *encrypted) {
    for (int i=0; i<7; i++) {
        pass[i] ^= 6-i;
        encrypted[i] = pass[i];
    }
}

void average_BGR(tagBITMAPINFOHEADER bitmapInfoHeader, byte *bitmap_BGR_data, int lower_bound, int upper_bound)
{
    for (int i = lower_bound; i < upper_bound; i += 3)
    {
        bitmap_BGR_data[i] = (bitmap_BGR_data[i] + bitmap_BGR_data[i + 1] + bitmap_BGR_data[i + 2]) / 3;
        bitmap_BGR_data[i + 1] = bitmap_BGR_data[i];
        bitmap_BGR_data[i + 2] = bitmap_BGR_data[i];
    }
}

void strconcat(char *str, char *parentcode, char add)
{
    int i = 0;
    while (*(parentcode + i) != '\0')
    {
        *(str + i) = *(parentcode + i);
        i++;
    }
    if (add != '2')
    {
        str[i] = add;
        str[i + 1] = '\0';
    }
    else
        str[i] = '\0';
}

void decimal_to_array_binary(int num, int *num_bin, int bits, int lower_bound) {
    int i=bits;
    while (num>=0 && i>lower_bound) {
        num_bin[i-1] = num%2;
        num /= 2;
        i--;
    }
}

void binary_to_literal_bits(int *num_bin, int bits, FILE *file) {
    int curr_bit = 0, bytes = bits/8;
    byte bit_buffer;
    if (bytes == 0) bytes = 1;

    for (int i=0; i<bits; i++) {
        if (num_bin[i])
            bit_buffer |= (1<<curr_bit);
        
        curr_bit++;
        if (curr_bit == 8) {
            fwrite(&bit_buffer, 1, 1, file);
            curr_bit = 0;
            bit_buffer = 0;
        }
    }
}

int binary_string_to_decimal(char *bin) {
    int val = 0;
    while (*bin != '\0')
        val = 2 * val + (*bin++ - '0');
    return val;
}

void literal_bits_to_binary_string(FILE *file, int bytes, char *str) {
    int lower_bound = 0, upper_bound = 8;
    for (int i=0; i<bytes; i++) {
        char temp_str[8];
        byte bits, temp, mask = 1;
        fread(&bits, 1, 1, file);
        for (int k=0; k<8; k++) {
            temp = bits & mask;
            if (temp == 1) {
                temp_str[k] = '1';
            }
            else {
                temp_str[k] = '0';
            }
            bits >>= 1;
        }
        int counter = 0;
        for (int j=lower_bound; j<upper_bound; j++) {
            if (temp_str[counter] == '1')
                str[j] = '1';
            else
                str[j] = '0';
            counter++;
        }
        lower_bound += 8;
        upper_bound += 8;
    }
}

void frequency_of_occurrence(tagBITMAPINFOHEADER info_header, int *frequency, unsigned char *bitmap_BGR_data)
{
    for (int i = 0; i < 256; i++)
        frequency[i] = 0;
    for (int i = 0; i < info_header.biSizeImage; i += 3)
        frequency[bitmap_BGR_data[i]] += 1;
}

int number_of_nodes(int *frequency)
{
    int nodes = 0;
    for (int i = 0; i < 256; i++)
        if (frequency[i] != 0)
            nodes++;
    return nodes;
}

void initialize_nodes(tagBITMAPINFOHEADER info_header, Node *pix_node, HuffTree *huff_node, int *frequency)
{
    int j = 0;
    float new_prob;
    for (int i = 0; i < 256; i++)
        if (frequency[i] != 0) {
            huff_node[j].pixel = i;
            pix_node[j].pixel = i;

            huff_node[j].arr_loc = j;

            new_prob = (float)frequency[i] / (float)info_header.biSizeImage;
            pix_node[j].freq = new_prob;
            huff_node[j].freq = new_prob;

            pix_node[j].left = NULL;
            pix_node[j].right = NULL;

            pix_node[j].huffman_code[0] = '\0';
            j++;
        }
}

void sort_node(Node *pix_node, HuffTree *huff_node, int nodes)
{
    // Sorting histogram
    HuffTree temp;
    for (int i = 0; i < nodes; i++)
        for (int j = i + 1; j < nodes; j++)
            if (huff_node[i].freq < huff_node[j].freq)
            {
                temp = huff_node[i];
                huff_node[i] = huff_node[j];
                huff_node[j] = temp;
            }
}

void build_huffman_tree(Node *pix_node, HuffTree *huff_node, int nodes)
{
    // Building Huffman Tree
    float sumprob;
    int sumpix;
    int i = 0, n = 0, k = 0;
    int next_node = nodes;

    while (n < nodes - 1) {
        // Adding the lowest two probabilities
        sumprob = huff_node[nodes - n - 1].freq + huff_node[nodes - n - 2].freq;
        sumpix = huff_node[nodes - n - 1].pixel + huff_node[nodes - n - 2].pixel;

        // Appending to the pix_node Array
        pix_node[next_node].pixel = sumpix;
        pix_node[next_node].freq = sumprob;
        pix_node[next_node].left = &pix_node[huff_node[nodes - n - 2].arr_loc];
        pix_node[next_node].right = &pix_node[huff_node[nodes - n - 1].arr_loc];
        pix_node[next_node].huffman_code[0] = '\0';
        i = 0;

        // keep sorting and update the huff_node array
        while (sumprob <= huff_node[i].freq)
            i++;

        //insert new node
        for (k = nodes; k >= 0; k--) {
            if (k == i) {
                huff_node[k].pixel = sumpix;
                huff_node[k].freq = sumprob;
                huff_node[k].arr_loc = next_node;
            }
            else if (k > i)
                huff_node[k] = huff_node[k - 1];
        }
        n += 1;
        next_node += 1;
    }
}

void assign_tree(Node *pix_node, int total_nodes, int nodes)
{
    char left = '0';
    char right = '1';
    for (int i = total_nodes - 1; i >= nodes; i--)
    {
        if (pix_node[i].left != NULL)
            strconcat(pix_node[i].left->huffman_code, pix_node[i].huffman_code, left);
        if (pix_node[i].right != NULL)
            strconcat(pix_node[i].right->huffman_code, pix_node[i].huffman_code, right);
    }
}

int total_huff_length(tagBITMAPINFOHEADER info_header, byte *bitmap_BGR_data, Node *pix_node, int nodes) {
    int val = 0, huff_length = 0;
    //Find total length of char array of all huffman code
    for (int i = 0; i < info_header.biSizeImage; i += 3) {
        val = (int)bitmap_BGR_data[i];
        for (int j = 0; j < nodes; j++)
            if (pix_node[j].pixel == val) {
                huff_length += strlen(pix_node[j].huffman_code);
            }
    }
    return huff_length;
}

void write_hbmp(tagBITMAPFILEHEADER file_header, tagBITMAPINFOHEADER info_header, char *file_name, byte *bitmap_BGR_data, int nodes, Node *pix_node, int *frequency, char *encrypted, char *flag)
{
    int val = 0, huff_length = total_huff_length(info_header, bitmap_BGR_data, pix_node, nodes), temp = 0, lower_bound = 0;
    int pass_bin[8] = {0}, nodes_bin[16] = {0}, val_bin[8] = {0}, freq_bin[24] = {0}, huff_length_bin[24] = {0}, flag_bin[8] = {0};
    
    if (strcmp(flag, "c") == 0) flag_bin[0] = 1;
    FILE *fptr = fopen(file_name, "wb");

    //Write File header and Info header
    fwrite(&file_header, sizeof(tagBITMAPFILEHEADER), 1, fptr);
    fwrite(&info_header, sizeof(tagBITMAPINFOHEADER), 1, fptr);

    //Write encrypted password
    fprintf(fptr, "%s", encrypted);

    //Write color flag (color(0x01) or grayscale(0x00))
    binary_to_literal_bits(flag_bin, 8, fptr);

    //Write total Huff Length (total bits + padding to make it multiple of 8)
    decimal_to_array_binary(huff_length, huff_length_bin, 24, 0); //WITHOUT Padding
    binary_to_literal_bits(huff_length_bin, 24, fptr);
    while (huff_length%8 != 0) //PADDING
        huff_length++; //make sure huff_length is a multiple of 8
    decimal_to_array_binary(huff_length, huff_length_bin, 24, 0); //WITH Padding
    binary_to_literal_bits(huff_length_bin, 24, fptr);
    
    //Write total leaf nodes in literal bits
    decimal_to_array_binary(nodes, nodes_bin, 16, 0);
    binary_to_literal_bits(nodes_bin, 16, fptr);

    //Write dictionary: value-freq in literal bits
    for (int i=0; i<256; i++) {
        if (frequency[i] != 0) {
            decimal_to_array_binary(i, val_bin, 8, 0);
            binary_to_literal_bits(val_bin, 8, fptr);
            decimal_to_array_binary(frequency[i], freq_bin, 24, 0);
            binary_to_literal_bits(freq_bin, 24, fptr);
        }
    }

    //Write huffman codes
    int *huff_bin = (int *)malloc(huff_length*sizeof(int));
    for (int i=0; i<huff_length; i++) huff_bin[i] = 0; //initialize huff_binary
    //Write huffman code in literal bits
    for (int i = 0; i < info_header.biSizeImage; i += 3) {
        val = (int)bitmap_BGR_data[i];
        for (int j = 0; j < nodes; j++)
            if (pix_node[j].pixel == val) {
                lower_bound = temp;
                temp += strlen(pix_node[j].huffman_code);
                decimal_to_array_binary(binary_string_to_decimal(pix_node[j].huffman_code), huff_bin, temp, lower_bound);
            }
    }
    //Write huffman code in a multiple of 8 literal bits
    int curr_bit = 0, upper_bound = 7;
    byte bit_buffer;
    lower_bound=0;
    for(int i=0; i<huff_length/8;i++) {
        for (int j=upper_bound; j>=lower_bound; j--) {
            if (huff_bin[j])
                bit_buffer |= (1<<curr_bit);
        
            curr_bit++;
            if (curr_bit == 8) {
                fwrite(&bit_buffer, 1, 1, fptr);
                curr_bit = 0;
                bit_buffer = 0;
            }
        } 
        upper_bound += 8;
        lower_bound +=8;
    }
    
    fclose(fptr);

    printf("'%s' file successfully created.\n", file_name);
}

FILE *load_hbmp_file(char *filename, tagBITMAPINFOHEADER *bitmapInfoHeader, tagBITMAPFILEHEADER *bitmapFileHeader, int *frequency, int *nodes, int *huff_length, int *huff_length_pure, char *encrypted, int *flag) {
    FILE *file = load_header(filename, bitmapInfoHeader, bitmapFileHeader);
    char *pass_bin, *nodes_bin, *val_bin, *freq_bin, *huff_length_bin, *huff_length_pure_bin;
    int val, freq;

    pass_bin = (char *)malloc(1*8*sizeof(char));
    nodes_bin = (char *)malloc(2*8*sizeof(char));
    val_bin = (char *)malloc(1*8*sizeof(char));
    huff_length_bin = (char *)malloc(3*8*sizeof(char));
    huff_length_pure_bin = (char *)malloc(3*8*sizeof(char));

    // move file pointer to point to dictionary (value, frequency)
    fseek(file, bitmapFileHeader->bfOffBits, SEEK_SET);

    //load encrypted password
    fread(encrypted, 7, 1, file);

    //load color flag (color(0x01) or grayscale(0x00))
    fread(flag, 1, 1, file);

    //load PURE Huff length (without padding)
    literal_bits_to_binary_string(file, 3, huff_length_pure_bin);
    *huff_length_pure = binary_string_to_decimal(huff_length_pure_bin);

    //load total Huff Length (total bits + padding to make it multiple of 8)
    literal_bits_to_binary_string(file, 3, huff_length_bin);
    *huff_length = binary_string_to_decimal(huff_length_bin);

    //read how many nodes (2 bytes)
    literal_bits_to_binary_string(file, 2, nodes_bin);
    *nodes = binary_string_to_decimal(nodes_bin);

    freq_bin = (char *)malloc(3*8*sizeof(char));
    //initialize frequency to 0
    for (int i = 0; i < 256; i++)
        frequency[i] = 0;
    //populate frequency
    for (int i=0; i<*nodes; i++) {
        literal_bits_to_binary_string(file, 1, val_bin);
        val = binary_string_to_decimal(val_bin);
        literal_bits_to_binary_string(file, 3, freq_bin);
        freq = binary_string_to_decimal(freq_bin);
        frequency[val] = freq;
    }
    return file;
}

void populate_BGR(FILE *file, tagBITMAPINFOHEADER info_header, byte *bitmap_BGR_data, Node *pix_node, int total_nodes, int nodes, int huff_length, int huff_length_pure, int flag) {
    int counter = 0, iter = 0; 
    char *huff_code = (char *)malloc(huff_length*sizeof(char));
    byte huff_byte, temp, mask = 1;
    // for (int i=0; i<huff_length; i++) huff_code[i] = '0'; //Initialize huff_code string

    //Store all huffman code into a string
    while(!feof(file)) {
        fread(&huff_byte, 1, 1, file);
        for (int i=7; i>=0; i--) {
            temp = huff_byte & mask;
            if (temp == 1)
                huff_code[counter+i] = '1';
            else
                huff_code[counter+i] = '0';
            huff_byte >>= 1;
        }
        counter += 8;
    }

    Node *root = &pix_node[total_nodes-1], *walker = root;

    if (flag == 1) {
        for (int i=0; i<huff_length; i++) {
            if (huff_code[i] == '0') {
                walker = walker->left;
                if (walker->left == NULL && walker->right == NULL) {
                    bitmap_BGR_data[iter] = (byte)walker->pixel;
                    walker = root;
                    iter++;
                }
            }
            if (huff_code[i] == '1') {
                walker = walker->right;
                if (walker->left == NULL && walker->right == NULL) {
                    bitmap_BGR_data[iter] = (byte)walker->pixel;
                    walker = root;
                    iter++;
                }
            }
        }
    }
    else {
        for (int i=0; i<huff_length_pure; i++) {
            if (huff_code[i] == '0') {
                walker = walker->left;
                if (walker->left == NULL && walker->right == NULL) {
                    bitmap_BGR_data[iter+0] = (byte)walker->pixel;
                    bitmap_BGR_data[iter+1] = (byte)walker->pixel;
                    bitmap_BGR_data[iter+2] = (byte)walker->pixel;
                    walker = root;
                    iter += 3;
                }
            }
            if (huff_code[i] == '1') {
                walker = walker->right;
                if (walker->left == NULL && walker->right == NULL) {
                    bitmap_BGR_data[iter+0] = (byte)walker->pixel;
                    bitmap_BGR_data[iter+1] = (byte)walker->pixel;
                    bitmap_BGR_data[iter+2] = (byte)walker->pixel;
                    walker = root;
                    iter += 3;
                }
            }
        }
    }
}

void write_bmp(tagBITMAPFILEHEADER file_header, tagBITMAPINFOHEADER info_header, byte *bitmap_BGR_data, char *file_name) {
    FILE *fptr = fopen(file_name, "wb");

    //Write File header and Info header
    fwrite(&file_header, sizeof(tagBITMAPFILEHEADER), 1, fptr);
    fwrite(&info_header, sizeof(tagBITMAPINFOHEADER), 1, fptr);
    fwrite(bitmap_BGR_data, info_header.biSizeImage, 1, fptr);

    fclose(fptr);
    printf("'%s' file successfully created.\n", file_name);
}

int main(int argc, char *argv[])
{
    int format, frequency[256], nodes, total_nodes;
    verify_input(argc, argv, &format);
    tagBITMAPFILEHEADER file_header;
    tagBITMAPINFOHEADER info_header;
    byte *bitmap_BGR_data;

    //COMPRESS
    if (format == 0) {
        clock_t start = clock();
        char pass[8], encrypted[8];
        strcpy(pass, argv[3]);

        printf("Compressing %s...\n", argv[2]);
        encrypt_and_decrypt(pass, encrypted);

        bitmap_BGR_data = load_bmp_file(argv[2], &info_header, &file_header);
        if (strcmp(argv[1], "g") == 0) {
            printf("Converting image to greyscale...\n");
            if (fork() == 0) {
                average_BGR(info_header, bitmap_BGR_data, 0, info_header.biSizeImage/2);
                return 0;
            }
            else {
                wait(0);
                average_BGR(info_header, bitmap_BGR_data, info_header.biSizeImage/2, info_header.biSizeImage);
            }
        }
        frequency_of_occurrence(info_header, frequency, bitmap_BGR_data);
        nodes = number_of_nodes(frequency);
        total_nodes = 2 * nodes - 1;
        Node *pix_node = (Node *)malloc(total_nodes * sizeof(Node));
        HuffTree *huff_node = (HuffTree *)malloc(nodes * sizeof(HuffTree));

        initialize_nodes(info_header, pix_node, huff_node, frequency);
        sort_node(pix_node, huff_node, nodes);
        build_huffman_tree(pix_node, huff_node, nodes);
        assign_tree(pix_node, total_nodes, nodes);

        write_hbmp(file_header, info_header, "out.hbmp", bitmap_BGR_data, nodes, pix_node, frequency, encrypted, argv[1]);
        clock_t end = clock();
        printf("Time: %f s\n", (double)(end-start)/CLOCKS_PER_SEC);
    }
    //DECOMPRESS
    else if (format == 1) {
        clock_t start = clock();
        FILE *file;
        int huff_length, huff_length_pure, flag = 0;
        char pass[8] = {0}, encrypted[8] = {0}, decrypted[8] = {0};
        strcpy(pass, argv[3]);
        
        file = load_hbmp_file(argv[2], &info_header, &file_header, frequency, &nodes, &huff_length, &huff_length_pure, decrypted, &flag);
        encrypt_and_decrypt(decrypted, encrypted);

        if (strcmp(pass, decrypted)) {
            printf("Password do not match with encrypted! Please try again.\n");
            exit(1);
        }
        printf("Decompressing %s...\n", argv[2]);

        bitmap_BGR_data = (byte *)malloc(info_header.biSizeImage);
        total_nodes = 2 * nodes - 1;
        Node *pix_node = (Node *)malloc(total_nodes * sizeof(Node));
        HuffTree *huff_node = (HuffTree *)malloc(nodes * sizeof(HuffTree));

        initialize_nodes(info_header, pix_node, huff_node, frequency);
        sort_node(pix_node, huff_node, nodes);
        build_huffman_tree(pix_node, huff_node, nodes);
        assign_tree(pix_node, total_nodes, nodes);

        populate_BGR(file, info_header, bitmap_BGR_data, pix_node, total_nodes, nodes, huff_length, huff_length_pure, flag);
        write_bmp(file_header, info_header, bitmap_BGR_data, "out.bmp");
        clock_t end = clock();
        printf("Time: %f s\n", (double)(end-start)/CLOCKS_PER_SEC);
    }
    
    return 0;
}
