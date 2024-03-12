
package operator.mutation;

import java.util.LinkedList;

/**
 *
 * @author al_x3
 */
public interface MutationOperator<T> {
    public T execute(T child);
}
