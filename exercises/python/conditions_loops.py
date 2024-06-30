if __name__ == "__main__":
    # Given: Two positive integers a and b (a<b<10000)
    # Return: The sum of all odd integers from a through b, inclusively.
    a, b = 4136, 8412
    sum = 0
    for i in range(a, b+1):
        if i % 2 != 0:
            sum += i
    print(sum)