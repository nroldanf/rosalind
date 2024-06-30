if __name__ == "__main__":
    input_string = "7nEKSikPnlDr1CapreolusuIdw6PbAdmlBC9SNkhppXVnp5YyHUcuikJsMzZP1fHkOeQ0lwD27WzlzTDbbNSPnigriceps653zT8r8CDIfY419TvRIpSKfwRAejMSknrsnzgapOzPjRqmfR8qJ0zfxT63snMMnc55iBSUF0SU7uo4"
    # Indices a, b, c, d
    indices = [13, 21, 85, 93]
    # Result: The slice of this string from indices a through b and c through d (inclusive)
    print(input_string[indices[0]:indices[1]+1], input_string[indices[2]:indices[3]+1])
