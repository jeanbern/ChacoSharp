namespace ChacoSharp.Utilities
{
    public static class BitReverse
    {
        /* Reverse the bits of a number. */
        public static int bit_reverse(int val, int nbits)
        /* value to reverse bits of */
        /* number of significant bits */
        {
            int mask_low, mask_high; /* masks for bits to interchange */
            int bit_low, bit_high; /* values of the bits in question */
            int i; /* loop counter */

            mask_low = 1;
            mask_high = 1 << (nbits - 1);
            for (i = 0; i < nbits / 2; i++)
            {
                bit_low = (val & mask_low) >> i;
                bit_high = (val & mask_high) >> (nbits - i - 1);
                mask_low <<= 1;
                mask_high >>= 1;
                if (bit_low != bit_high)
                {
                    val ^= (1 << i) ^ (1 << (nbits - i - 1));
                }
            }

            return (val);
        }
    }
}
