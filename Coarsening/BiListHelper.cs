namespace ChacoSharp.Coarsening
{
    public static unsafe class BiListHelper
    {
        /* Note: bi-directional lists aren't assumed to be sorted. */

        public static void add2bilist( /* add val to unsorted list */
            bilist* lptr, /* element to add */
            bilist** list /* list added to */
        )
        {
            lptr->next = *list;
            if (*list != null)
            {
                (*list)->prev = lptr;
            }

            lptr->prev = null;
            *list = lptr;
        }

        public static void removebilist(bilist* lptr, /* ptr to element to remove */
            bilist** list /* head of list to remove it from */
        )

/* Remove an element from a bidirectional list. */
        {
            if (lptr->next != null)
            {
                lptr->next->prev = lptr->prev;
            }

            if (lptr->prev != null)
            {
                lptr->prev->next = lptr->next;
            }
            else
            {
                *list = lptr->next;
            }
        }

        public static void movebilist(bilist* lptr, /* ptr to element to move */
            bilist** oldlist, /* head of list to remove it from */
            bilist** newlist /* head of list to add it to */
        )

/* Move an element from a old bidirectional list to new one. */
        {
            removebilist(lptr, oldlist);

            add2bilist(lptr, newlist);
        }
    }
}
